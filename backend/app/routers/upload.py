"""Upload router — accepts count matrix and optional metadata files."""
from __future__ import annotations

import json
import uuid
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, File, HTTPException, UploadFile

router = APIRouter()

_MAX_SIZE = 50 * 1024 * 1024  # 50 MB
_ALLOWED_EXTS = {".tsv", ".csv", ".txt", ".xlsx"}
_UPLOAD_DIR = Path("/tmp/bioflowvalidator/uploads")

# Internal storage names — never user-provided, preventing path traversal
_COUNT_STORAGE_NAME = "count.file"
_META_STORAGE_NAME = "meta.file"
_MANIFEST_NAME = "manifest.json"


def _validate_upload(file: UploadFile) -> None:
    ext = Path(file.filename or "").suffix.lower()
    if ext not in _ALLOWED_EXTS:
        raise HTTPException(
            status_code=400,
            detail=f"File '{file.filename}' has unsupported extension '{ext}'. "
            f"Allowed: {sorted(_ALLOWED_EXTS)}",
        )


@router.post("")
async def upload_files(
    count_matrix: UploadFile = File(..., description="Count matrix file (TSV/CSV)"),
    metadata: Optional[UploadFile] = File(None, description="Sample metadata file (TSV/CSV)"),
) -> dict:
    """Accept uploaded files, persist them under fixed internal names, return a job_id."""
    _validate_upload(count_matrix)

    job_id = str(uuid.uuid4())
    job_dir = _UPLOAD_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    # Register the trusted path immediately — validate router will look it up
    from app import store as _store
    _store.register_job_dir(job_id, job_dir)

    chunk_size = 1024 * 1024  # 1MB
    try:
        # Read count_matrix in chunks to check size
        count_data = b""
        while True:
            chunk = await count_matrix.read(chunk_size)
            if not chunk:
                break
            count_data += chunk
            if len(count_data) > _MAX_SIZE:
                raise HTTPException(status_code=413, detail="Count matrix file exceeds 50 MB limit.")

        # Store with a fixed, safe internal name
        (job_dir / _COUNT_STORAGE_NAME).write_bytes(count_data)

        manifest: dict = {
            "count_matrix_filename": count_matrix.filename or "counts.tsv",
            "count_matrix_size": len(count_data),
        }

        if metadata and metadata.filename:
            _validate_upload(metadata)
            # Read metadata in chunks to check size
            meta_data = b""
            while True:
                chunk = await metadata.read(chunk_size)
                if not chunk:
                    break
                meta_data += chunk
                if len(meta_data) > _MAX_SIZE:
                    raise HTTPException(status_code=413, detail="Metadata file exceeds 50 MB limit.")

            (job_dir / _META_STORAGE_NAME).write_bytes(meta_data)
            manifest["metadata_filename"] = metadata.filename
            manifest["metadata_size"] = len(meta_data)

        (job_dir / _MANIFEST_NAME).write_text(json.dumps(manifest))

    except Exception:
        # Clean up directory on failure
        if job_dir.exists():
            import shutil
            shutil.rmtree(job_dir)
        raise

    return {"job_id": job_id, **manifest}
