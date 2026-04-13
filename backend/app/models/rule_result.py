"""Data models for rule results and validation reports."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Literal


RuleSeverity = Literal["ERROR", "WARNING", "INFO"]
RuleStatus = Literal["PASS", "FAIL", "SKIP"]


@dataclass
class RuleResult:
    rule_id: str
    category: str
    severity: RuleSeverity
    status: RuleStatus
    message: str
    affected_items: list[str] = field(default_factory=list)
    suggestion: str = ""
    details: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "rule_id": self.rule_id,
            "category": self.category,
            "severity": self.severity,
            "status": self.status,
            "message": self.message,
            "affected_items": self.affected_items,
            "suggestion": self.suggestion,
            "details": self.details,
        }


@dataclass
class FileMeta:
    filename: str
    size_bytes: int
    sha256: str

    @staticmethod
    def from_bytes(filename: str, data: bytes) -> "FileMeta":
        return FileMeta(
            filename=filename,
            size_bytes=len(data),
            sha256=hashlib.sha256(data).hexdigest(),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "filename": self.filename,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
        }


@dataclass
class ValidationReport:
    job_id: str
    timestamp: str
    files: list[FileMeta]
    results: list[RuleResult]
    error_count: int = 0
    warning_count: int = 0
    pass_count: int = 0
    skip_count: int = 0

    def __post_init__(self) -> None:
        self._recount()

    def _recount(self) -> None:
        self.error_count = sum(
            1 for r in self.results if r.status == "FAIL" and r.severity == "ERROR"
        )
        self.warning_count = sum(
            1 for r in self.results if r.status == "FAIL" and r.severity == "WARNING"
        )
        self.pass_count = sum(1 for r in self.results if r.status == "PASS")
        self.skip_count = sum(1 for r in self.results if r.status == "SKIP")

    @staticmethod
    def build(
        job_id: str, files: list[FileMeta], results: list[RuleResult]
    ) -> "ValidationReport":
        ts = datetime.now(timezone.utc).isoformat()
        report = ValidationReport(
            job_id=job_id, timestamp=ts, files=files, results=results
        )
        return report

    def to_dict(self) -> dict[str, Any]:
        return {
            "job_id": self.job_id,
            "timestamp": self.timestamp,
            "files": [f.to_dict() for f in self.files],
            "summary": {
                "error_count": self.error_count,
                "warning_count": self.warning_count,
                "pass_count": self.pass_count,
                "skip_count": self.skip_count,
                "total_rules": len(self.results),
            },
            "results": [r.to_dict() for r in self.results],
        }
