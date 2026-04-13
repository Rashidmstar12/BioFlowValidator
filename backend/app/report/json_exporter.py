"""JSONExporter — serialises a ValidationReport to JSON."""
from __future__ import annotations

import json

from app.models.rule_result import ValidationReport


def to_json(report: ValidationReport, indent: int = 2) -> str:
    return json.dumps(report.to_dict(), indent=indent, ensure_ascii=False)
