"""Validate router — triggers the validation engine for a given job."""
from __future__ import annotations

import json
import uuid
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, field_validator

from app import store
from app.engine.parser import parse_files
from app.engine.runner import run_all
from app.models.rule_result import FileMeta
from app.report.builder import build_report

router = APIRouter()

_UPLOAD_DIR = Path("/tmp/bioflowvalidator/uploads")

# These are the fixed internal storage names written by the upload router.
# They are never derived from user input, preventing path traversal.
_COUNT_STORAGE_NAME = "count.file"
_META_STORAGE_NAME = "meta.file"
_MANIFEST_NAME = "manifest.json"


class ValidateRequest(BaseModel):
    job_id: str

    @field_validator("job_id")
    @classmethod
    def _validate_job_id(cls, v: str) -> str:
        """Reject non-UUID job IDs to prevent path traversal."""
        try:
            uuid.UUID(v)
        except ValueError:
            raise ValueError("job_id must be a valid UUID")
        return v


@router.post("")
async def validate(req: ValidateRequest) -> dict:
    """Run validation engine on previously uploaded files."""
    # Only UUID job_ids reach here; safe to use directly in path construction.
    job_dir = _UPLOAD_DIR / req.job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail=f"Job '{req.job_id}' not found. Upload files first.")

    manifest_path = job_dir / _MANIFEST_NAME
    if not manifest_path.exists():
        raise HTTPException(status_code=404, detail="Job manifest not found.")
    manifest = json.loads(manifest_path.read_text())

    count_path = job_dir / _COUNT_STORAGE_NAME
    if not count_path.exists():
        raise HTTPException(status_code=404, detail="Count matrix file not found in job directory.")

    count_filename: str = manifest.get("count_matrix_filename", "counts.tsv")
    count_bytes = count_path.read_bytes()
    file_metas = [FileMeta.from_bytes(count_filename, count_bytes)]

    meta_bytes: Optional[bytes] = None
    meta_filename = ""
    meta_path = job_dir / _META_STORAGE_NAME
    if meta_path.exists():
        meta_bytes = meta_path.read_bytes()
        meta_filename = manifest.get("metadata_filename", "metadata.tsv")
        file_metas.append(FileMeta.from_bytes(meta_filename, meta_bytes))

    context = parse_files(
        count_bytes=count_bytes,
        count_filename=count_filename,
        metadata_bytes=meta_bytes,
        metadata_filename=meta_filename,
    )

    results = run_all(context)
    report = build_report(
        job_id=req.job_id,
        files=file_metas,
        results=results,
    )
    store.save(req.job_id, report)

    return {
        "job_id": req.job_id,
        "status": "complete",
        "summary": report.to_dict()["summary"],
    }
