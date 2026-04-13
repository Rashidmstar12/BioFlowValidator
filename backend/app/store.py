"""Shared in-memory job store (MVP — replace with Redis/DB for production)."""
from __future__ import annotations

from app.models.rule_result import ValidationReport

_jobs: dict[str, ValidationReport] = {}


def save(job_id: str, report: ValidationReport) -> None:
    _jobs[job_id] = report


def get(job_id: str) -> ValidationReport | None:
    return _jobs.get(job_id)
