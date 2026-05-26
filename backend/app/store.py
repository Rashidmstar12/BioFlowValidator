"""Shared persistent store for both validation reports and upload job directories."""
from __future__ import annotations

import json
import os
from pathlib import Path

from app.models.rule_result import FileMeta, RuleResult, ValidationReport

_UPLOAD_DIR = Path("/tmp/bioflowvalidator/uploads")
_REGISTRY_FILE = _UPLOAD_DIR / "registry.json"

# In-memory caches
_reports: dict[str, ValidationReport] = {}
_job_dirs: dict[str, Path] = {}


def _load_registry() -> dict[str, str]:
    if not _REGISTRY_FILE.exists():
        return {}
    try:
        with open(_REGISTRY_FILE, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def _save_registry(registry: dict[str, str]) -> None:
    _UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    try:
        temp_file = _REGISTRY_FILE.with_suffix(".tmp")
        with open(temp_file, "w", encoding="utf-8") as f:
            json.dump(registry, f, indent=2)
        if temp_file.exists():
            if _REGISTRY_FILE.exists():
                os.remove(_REGISTRY_FILE)
            os.rename(temp_file, _REGISTRY_FILE)
    except Exception:
        pass


# ── Report store ──────────────────────────────────────────────────────────────

def save(job_id: str, report: ValidationReport) -> None:
    # Save to memory cache
    _reports[job_id] = report

    # Save to disk inside job_dir as report.json
    job_dir = get_job_dir(job_id)
    if job_dir:
        job_dir.mkdir(parents=True, exist_ok=True)
        report_file = job_dir / "report.json"
        try:
            with open(report_file, "w", encoding="utf-8") as f:
                json.dump(report.to_dict(), f, indent=2)
        except Exception:
            pass


def get(job_id: str) -> ValidationReport | None:
    # Check memory cache first
    if job_id in _reports:
        return _reports[job_id]

    # Look up job_dir from registry
    job_dir = get_job_dir(job_id)
    if job_dir:
        report_file = job_dir / "report.json"
        if report_file.exists():
            try:
                with open(report_file, "r", encoding="utf-8") as f:
                    data = json.load(f)

                # Reconstruct ValidationReport object
                files = [FileMeta(f["filename"], f["size_bytes"], f["sha256"]) for f in data["files"]]
                results = [
                    RuleResult(
                        rule_id=r["rule_id"],
                        category=r["category"],
                        severity=r["severity"],
                        status=r["status"],
                        message=r["message"],
                        affected_items=r.get("affected_items", []),
                        suggestion=r.get("suggestion", ""),
                        details=r.get("details", {})
                    )
                    for r in data["results"]
                ]
                report = ValidationReport(
                    job_id=data["job_id"],
                    timestamp=data["timestamp"],
                    files=files,
                    results=results,
                    error_count=data["summary"]["error_count"],
                    warning_count=data["summary"]["warning_count"],
                    pass_count=data["summary"]["pass_count"],
                    skip_count=data["summary"]["skip_count"]
                )
                # Cache it
                _reports[job_id] = report
                return report
            except Exception:
                return None
    return None


# ── Job directory registry ────────────────────────────────────────────────────

def register_job_dir(job_id: str, job_dir: Path) -> None:
    """Called by the upload router to register a trusted path for a job."""
    _job_dirs[job_id] = job_dir

    registry = _load_registry()
    registry[job_id] = str(job_dir.resolve())
    _save_registry(registry)


def get_job_dir(job_id: str) -> Path | None:
    """Return the trusted job directory for a job_id, or None if unknown."""
    if job_id in _job_dirs:
        return _job_dirs[job_id]

    registry = _load_registry()
    path_str = registry.get(job_id)
    if path_str:
        path = Path(path_str)
        if path.exists():
            _job_dirs[job_id] = path
            return path
    return None
