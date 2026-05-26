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


def test_full_validation_valid_excel():
    import pandas as pd
    import io

    # Build valid dataframes
    count_df = pd.DataFrame(
        {
            "ctrl_1": [10000, 5000, 30000],
            "ctrl_2": [11000, 5500, 29000],
            "ctrl_3": [9000, 4500, 31000],
            "treat_1": [20000, 15000, 10000],
            "treat_2": [21000, 14500, 10500],
            "treat_3": [19000, 15500, 9500]
        },
        index=["ENSG00000000001", "ENSG00000000002", "ENSG00000000003"]
    )
    meta_df = pd.DataFrame(
        {
            "condition": ["control", "control", "control", "treated", "treated", "treated"]
        },
        index=["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"]
    )
    
    count_out = io.BytesIO()
    count_df.to_excel(count_out, index=True)
    count_bytes = count_out.getvalue()

    meta_out = io.BytesIO()
    meta_df.to_excel(meta_out, index=True)
    meta_bytes = meta_out.getvalue()

    # Upload Excel files
    upload_r = client.post(
        "/upload",
        files={
            "count_matrix": ("counts.xlsx", io.BytesIO(count_bytes), "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
            "metadata": ("metadata.xlsx", io.BytesIO(meta_bytes), "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
        },
    )
    assert upload_r.status_code == 200
    job_id = upload_r.json()["job_id"]

    # Validate
    val_r = client.post("/validate", json={"job_id": job_id})
    assert val_r.status_code == 200
    summary = val_r.json()["summary"]
    # Excel files should parse and run correctly without errors on valid data
    assert summary["error_count"] == 0

    # Confirm html report includes details
    html_r = client.get(f"/report/{job_id}")
    assert html_r.status_code == 200
    assert "Show Rule Details" in html_r.text
