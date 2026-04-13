"""Validate router — triggers the validation engine for a given job."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from app import store
from app.engine.parser import parse_files
from app.engine.runner import run_all
from app.models.rule_result import FileMeta
from app.report.builder import build_report

router = APIRouter()

_UPLOAD_DIR = Path("/tmp/bioflowvalidator/uploads")


class ValidateRequest(BaseModel):
    job_id: str
    count_matrix_filename: str
    metadata_filename: Optional[str] = None


@router.post("")
async def validate(req: ValidateRequest) -> dict:
    """Run validation engine on previously uploaded files."""
    job_dir = _UPLOAD_DIR / req.job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail=f"Job '{req.job_id}' not found. Upload files first.")

    count_path = job_dir / req.count_matrix_filename
    if not count_path.exists():
        raise HTTPException(status_code=404, detail=f"Count matrix file not found in job directory.")

    count_bytes = count_path.read_bytes()
    file_metas = [FileMeta.from_bytes(req.count_matrix_filename, count_bytes)]

    meta_bytes: Optional[bytes] = None
    meta_filename = ""
    if req.metadata_filename:
        meta_path = job_dir / req.metadata_filename
        if meta_path.exists():
            meta_bytes = meta_path.read_bytes()
            file_metas.append(FileMeta.from_bytes(req.metadata_filename, meta_bytes))
            meta_filename = req.metadata_filename

    context = parse_files(
        count_bytes=count_bytes,
        count_filename=req.count_matrix_filename,
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
