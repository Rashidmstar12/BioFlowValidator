"""ReportBuilder — assembles a ValidationReport from rule results."""
from __future__ import annotations

from app.models.rule_result import FileMeta, RuleResult, ValidationReport


def build_report(
    job_id: str,
    files: list[FileMeta],
    results: list[RuleResult],
) -> ValidationReport:
    return ValidationReport.build(job_id=job_id, files=files, results=results)
