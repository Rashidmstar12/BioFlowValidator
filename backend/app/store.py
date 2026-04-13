"""Shared in-memory store for both validation reports and upload job directories."""
from __future__ import annotations

from pathlib import Path

from app.models.rule_result import ValidationReport

# Maps job_id → ValidationReport
_reports: dict[str, ValidationReport] = {}

# Maps job_id → trusted job directory Path object (set by upload router at upload time)
_job_dirs: dict[str, Path] = {}


# ── Report store ──────────────────────────────────────────────────────────────

def save(job_id: str, report: ValidationReport) -> None:
    _reports[job_id] = report


def get(job_id: str) -> ValidationReport | None:
    return _reports.get(job_id)


# ── Job directory registry ────────────────────────────────────────────────────

def register_job_dir(job_id: str, job_dir: Path) -> None:
    """Called by the upload router to register a trusted path for a job."""
    _job_dirs[job_id] = job_dir


def get_job_dir(job_id: str) -> Path | None:
    """Return the trusted job directory for a job_id, or None if unknown."""
    return _job_dirs.get(job_id)
