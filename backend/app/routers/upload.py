"""Upload router — accepts count matrix and optional metadata files."""
from __future__ import annotations

import uuid
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, File, HTTPException, UploadFile

router = APIRouter()

_MAX_SIZE = 50 * 1024 * 1024  # 50 MB
_ALLOWED_EXTS = {".tsv", ".csv", ".txt", ".xlsx"}
_UPLOAD_DIR = Path("/tmp/bioflowvalidator/uploads")


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
    """Accept uploaded files, persist them temporarily, return a job_id."""
    _validate_upload(count_matrix)

    job_id = str(uuid.uuid4())
    job_dir = _UPLOAD_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    count_data = await count_matrix.read()
    if len(count_data) > _MAX_SIZE:
        raise HTTPException(status_code=413, detail="Count matrix file exceeds 50 MB limit.")
    (job_dir / (count_matrix.filename or "counts.tsv")).write_bytes(count_data)

    meta_info: dict = {}
    if metadata and metadata.filename:
        _validate_upload(metadata)
        meta_data = await metadata.read()
        if len(meta_data) > _MAX_SIZE:
            raise HTTPException(status_code=413, detail="Metadata file exceeds 50 MB limit.")
        (job_dir / metadata.filename).write_bytes(meta_data)
        meta_info = {"metadata_filename": metadata.filename, "metadata_size": len(meta_data)}

    return {
        "job_id": job_id,
        "count_matrix_filename": count_matrix.filename,
        "count_matrix_size": len(count_data),
        **meta_info,
    }
