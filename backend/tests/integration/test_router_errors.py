"""Integration tests for router error paths not covered by the main integration suite."""
from __future__ import annotations

import io
import uuid

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

VALID_COUNTS = (
    "gene_id\tctrl_1\tctrl_2\tctrl_3\ttreat_1\ttreat_2\ttreat_3\n"
    "ENSG00000000001\t10000\t11000\t9000\t20000\t21000\t19000\n"
    "ENSG00000000002\t5000\t5500\t4500\t15000\t14500\t15500\n"
)

VALID_META = (
    "sample_id\tcondition\n"
    "ctrl_1\tcontrol\n"
    "ctrl_2\tcontrol\n"
    "ctrl_3\tcontrol\n"
    "treat_1\ttreated\n"
    "treat_2\ttreated\n"
    "treat_3\ttreated\n"
)


# ---------------------------------------------------------------------------
# Upload — oversized count matrix (line 51 in upload.py)
# ---------------------------------------------------------------------------

def test_upload_oversized_count_matrix(monkeypatch):
    """Files exceeding 50 MB must be rejected with 413."""
    from app.routers import upload as _upload

    original_max = _upload._MAX_SIZE
    monkeypatch.setattr(_upload, "_MAX_SIZE", 5)  # set limit to 5 bytes

    r = client.post(
        "/upload",
        files={
            "count_matrix": (
                "counts.tsv",
                io.BytesIO(VALID_COUNTS.encode()),
                "text/plain",
            )
        },
    )
    assert r.status_code == 413
    assert "50 MB" in r.json()["detail"] or "limit" in r.json()["detail"].lower()

    monkeypatch.setattr(_upload, "_MAX_SIZE", original_max)


# ---------------------------------------------------------------------------
# Upload — oversized metadata file (line 65 in upload.py)
# ---------------------------------------------------------------------------

def test_upload_oversized_metadata(monkeypatch):
    """Metadata file exceeding 50 MB must be rejected with 413."""
    from app.routers import upload as _upload

    # Set limit small enough that metadata bytes exceed it
    monkeypatch.setattr(_upload, "_MAX_SIZE", len(VALID_COUNTS.encode()) + 100)

    r = client.post(
        "/upload",
        files={
            "count_matrix": (
                "counts.tsv",
                io.BytesIO(VALID_COUNTS.encode()),
                "text/plain",
            ),
            "metadata": (
                "metadata.tsv",
                io.BytesIO(VALID_META.encode() * 100),
                "text/plain",
            ),
        },
    )
    assert r.status_code == 413


# ---------------------------------------------------------------------------
# Report router — HTML 404 (line 28 in report.py)
# ---------------------------------------------------------------------------

def test_html_report_unknown_job():
    """GET /report/<unknown_job> must return 404."""
    unknown_id = str(uuid.uuid4())
    r = client.get(f"/report/{unknown_id}")
    assert r.status_code == 404


# ---------------------------------------------------------------------------
# Validate router — 404 paths (lines 50, 57, 62 in validate.py)
# ---------------------------------------------------------------------------

def test_validate_valid_uuid_not_uploaded():
    """A valid UUID4 that was never uploaded must return 404."""
    fresh_id = str(uuid.uuid4())
    r = client.post("/validate", json={"job_id": fresh_id})
    assert r.status_code == 404


def test_validate_missing_manifest(tmp_path, monkeypatch):
    """Job dir exists but manifest.json is absent → 404."""
    from app import store as _store

    job_id = str(uuid.uuid4())
    job_dir = tmp_path / job_id
    job_dir.mkdir()
    _store.register_job_dir(job_id, job_dir)
    # Do NOT write manifest.json

    r = client.post("/validate", json={"job_id": job_id})
    assert r.status_code == 404


def test_validate_missing_count_file(tmp_path, monkeypatch):
    """Manifest exists but count.file is absent → 404."""
    import json as _json
    from app import store as _store

    job_id = str(uuid.uuid4())
    job_dir = tmp_path / job_id
    job_dir.mkdir()
    _store.register_job_dir(job_id, job_dir)

    manifest = {"count_matrix_filename": "counts.tsv", "count_matrix_size": 100}
    (job_dir / "manifest.json").write_text(_json.dumps(manifest))
    # Do NOT write count.file

    r = client.post("/validate", json={"job_id": job_id})
    assert r.status_code == 404
