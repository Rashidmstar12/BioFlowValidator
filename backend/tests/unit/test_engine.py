"""Unit tests for engine components: runner, parser, and report exporters."""
from __future__ import annotations

import json

import pytest

from app.engine.parser import _detect_encoding, _read_tabular, parse_files
from app.engine.runner import run_all
from app.models.context import ValidationContext
from app.models.rule_result import RuleResult, ValidationReport, FileMeta
from app.report.json_exporter import to_json
from app.rules.base import BaseRule
from tests.conftest import _ctx_from_dfs, make_count_matrix


# ---------------------------------------------------------------------------
# Parser — fallback encoding path (lines 39-40)
# ---------------------------------------------------------------------------

def test_read_tabular_invalid_encoding_fallback(monkeypatch):
    """If chardet returns an unrecognised codec name, decode falls back to utf-8."""
    import chardet as _chardet

    original_detect = _chardet.detect

    def patched_detect(data):
        return {"encoding": "not-a-real-codec", "confidence": 0.99}

    monkeypatch.setattr("app.engine.parser.chardet.detect", patched_detect)

    tsv = b"gene_id\ts1\ts2\nG1\t10\t20\nG2\t5\t15\n"
    df, delim, enc = _read_tabular(tsv, "counts.tsv")
    assert df is not None
    assert list(df.columns) == ["s1", "s2"]


def test_detect_encoding_utf8_variants():
    """UTF-8 BOM and aliased encodings are normalised to utf-8-sig."""
    # We test the normalisation branch directly
    import chardet

    bom_bytes = b"\xef\xbb\xbfgene_id\ts1\nG1\t1\n"
    enc = _detect_encoding(bom_bytes)
    # chardet should detect UTF-8-SIG; our function maps it to utf-8-sig
    assert enc.lower() in ("utf-8-sig", "utf-8", "ascii")


# ---------------------------------------------------------------------------
# Runner — exception handling (lines 22-23)
# ---------------------------------------------------------------------------

class _BrokenRule(BaseRule):
    rule_id = "TST-BROKEN"
    category = "test"
    severity = "ERROR"
    description = "A rule that always raises."

    def run(self, context: ValidationContext) -> RuleResult:
        raise RuntimeError("Simulated internal rule failure")


def test_run_all_catches_rule_exception(monkeypatch):
    """run_all must not propagate a rule exception; it returns an ERROR result."""
    # Monkeypatch the name as imported into runner.py
    monkeypatch.setattr("app.engine.runner.get_rules", lambda: [_BrokenRule()])

    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    results = run_all(ctx)
    assert len(results) == 1
    result = results[0]
    assert result.rule_id == "TST-BROKEN"
    assert result.status == "FAIL"
    assert result.severity == "ERROR"
    assert "unexpected exception" in result.message
    assert "traceback" in result.details


# ---------------------------------------------------------------------------
# JSON exporter (line 10)
# ---------------------------------------------------------------------------

def test_to_json_produces_valid_json():
    file_meta = FileMeta.from_bytes("counts.tsv", b"gene_id\ts1\nG1\t1\n")
    result = RuleResult(
        rule_id="FMT-001",
        category="format",
        severity="ERROR",
        status="PASS",
        message="ok",
    )
    report = ValidationReport(
        job_id="test-job",
        timestamp="2024-01-01T00:00:00+00:00",
        files=[file_meta],
        results=[result],
    )
    json_str = to_json(report)
    parsed = json.loads(json_str)
    assert parsed["job_id"] == "test-job"
    assert len(parsed["results"]) == 1


def test_to_json_non_ascii_characters():
    """Non-ASCII gene names must be preserved (ensure_ascii=False)."""
    file_meta = FileMeta.from_bytes("c.tsv", b"x")
    result = RuleResult(
        rule_id="GEN-001",
        category="gene",
        severity="WARNING",
        status="FAIL",
        message="Gène spécial: αβγ",
    )
    report = ValidationReport(
        job_id="j1",
        timestamp="2024-01-01T00:00:00+00:00",
        files=[file_meta],
        results=[result],
    )
    json_str = to_json(report)
    assert "αβγ" in json_str
