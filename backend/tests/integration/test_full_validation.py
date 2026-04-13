"""Integration tests — full API workflow."""
from __future__ import annotations

import io
import json

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def _tsv(df_str: str) -> bytes:
    return df_str.encode()


VALID_COUNTS = (
    "gene_id\tctrl_1\tctrl_2\tctrl_3\ttreat_1\ttreat_2\ttreat_3\n"
    "ENSG00000000001\t10000\t11000\t9000\t20000\t21000\t19000\n"
    "ENSG00000000002\t5000\t5500\t4500\t15000\t14500\t15500\n"
    "ENSG00000000003\t30000\t29000\t31000\t10000\t10500\t9500\n"
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


def test_health():
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["status"] == "ok"


def test_upload_count_only():
    r = client.post(
        "/upload",
        files={"count_matrix": ("counts.tsv", io.BytesIO(VALID_COUNTS.encode()), "text/plain")},
    )
    assert r.status_code == 200
    data = r.json()
    assert "job_id" in data
    assert data["count_matrix_filename"] == "counts.tsv"


def test_upload_with_metadata():
    r = client.post(
        "/upload",
        files={
            "count_matrix": ("counts.tsv", io.BytesIO(VALID_COUNTS.encode()), "text/plain"),
            "metadata": ("metadata.tsv", io.BytesIO(VALID_META.encode()), "text/plain"),
        },
    )
    assert r.status_code == 200
    assert "metadata_filename" in r.json()


def test_upload_rejects_bad_extension():
    r = client.post(
        "/upload",
        files={"count_matrix": ("counts.exe", io.BytesIO(b"bad"), "application/octet-stream")},
    )
    assert r.status_code == 400


def test_full_validation_valid_data():
    # Upload
    upload_r = client.post(
        "/upload",
        files={
            "count_matrix": ("counts.tsv", io.BytesIO(VALID_COUNTS.encode()), "text/plain"),
            "metadata": ("metadata.tsv", io.BytesIO(VALID_META.encode()), "text/plain"),
        },
    )
    job_id = upload_r.json()["job_id"]

    # Validate
    val_r = client.post("/validate", json={"job_id": job_id})
    assert val_r.status_code == 200
    summary = val_r.json()["summary"]
    assert summary["error_count"] == 0

    # Get JSON results
    res_r = client.get(f"/report/results/{job_id}")
    assert res_r.status_code == 200
    report = res_r.json()
    assert report["job_id"] == job_id
    assert len(report["results"]) > 0

    # Get HTML report
    html_r = client.get(f"/report/{job_id}")
    assert html_r.status_code == 200
    assert "BioFlowValidator" in html_r.text


def test_full_validation_sample_mismatch():
    mismatch_meta = (
        "sample_id\tcondition\n"
        "ctrl_1\tcontrol\n"
        "ctrl_2\tcontrol\n"
        "ctrl_3\tcontrol\n"
        "treat_1\ttreated\n"
        "treat_2\ttreated\n"
        # treat_3 missing
    )
    upload_r = client.post(
        "/upload",
        files={
            "count_matrix": ("counts.tsv", io.BytesIO(VALID_COUNTS.encode()), "text/plain"),
            "metadata": ("metadata.tsv", io.BytesIO(mismatch_meta.encode()), "text/plain"),
        },
    )
    job_id = upload_r.json()["job_id"]
    val_r = client.post("/validate", json={"job_id": job_id})
    assert val_r.json()["summary"]["error_count"] >= 1


def test_validate_unknown_job():
    r = client.post("/validate", json={"job_id": "does-not-exist"})
    # Invalid UUID → 422 (validation error) or 404; both are acceptable rejections
    assert r.status_code in (404, 422)


def test_results_unknown_job():
    r = client.get("/report/results/does-not-exist")
    assert r.status_code == 404
