"""Format validation rules."""
from __future__ import annotations

import re

import chardet
import numpy as np
import pandas as pd

from app.models.context import ValidationContext
from app.models.rule_result import RuleResult
from app.rules.base import BaseRule


class EncodingRule(BaseRule):
    rule_id = "FMT-001"
    category = "format"
    severity = "ERROR"
    description = "Count matrix file must be UTF-8 encoded."

    def run(self, context: ValidationContext) -> RuleResult:
        data = context.count_matrix_bytes
        if not data:
            return self._skip("No count matrix bytes available.")
        result = chardet.detect(data)
        enc = (result.get("encoding") or "").lower()
        if enc in ("", "ascii", "utf-8", "utf-8-sig", "utf8", "utf8sig", "utf-8sig"):
            return self._pass(f"File encoding detected as '{enc}' (UTF-8 compatible).")
        return self._fail(
            message=f"File encoding detected as '{enc}'. UTF-8 is required.",
            affected_items=[context.count_filename],
            suggestion=(
                "Re-save the file with UTF-8 encoding. "
                "In Excel: File → Save As → CSV UTF-8."
            ),
        )


class DelimiterRule(BaseRule):
    rule_id = "FMT-002"
    category = "format"
    severity = "WARNING"
    description = "Count matrix should use a consistent delimiter (TSV or CSV)."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix could not be parsed.")
        delim = context.count_delimiter
        name = "tab-separated (TSV)" if delim == "\t" else "comma-separated (CSV)"
        return self._pass(f"Delimiter detected as {name}.")


class HeaderRule(BaseRule):
    rule_id = "FMT-003"
    category = "format"
    severity = "ERROR"
    description = "Count matrix must have a header row with sample IDs."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix could not be parsed.")
        df = context.count_matrix
        if df.shape[1] == 0:
            return self._fail(
                "Count matrix has no columns. Header row may be missing.",
                suggestion="Ensure the first row contains sample IDs.",
            )
        # Heuristic: if all column names look like integers, the header is probably missing
        if all(str(c).isdigit() for c in df.columns):
            return self._fail(
                "All column names are numeric integers, which suggests the header row is missing.",
                affected_items=[str(c) for c in df.columns[:5]],
                suggestion="Add a header row with sample IDs as the first row.",
            )
        return self._pass("Header row with sample IDs detected.")


class DuplicateColumnRule(BaseRule):
    rule_id = "FMT-004"
    category = "format"
    severity = "ERROR"
    description = "Count matrix must not contain duplicate sample (column) names."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix could not be parsed.")
        cols = list(context.count_matrix.columns)
        seen: dict = {}
        dupes = []
        for c in cols:
            seen[c] = seen.get(c, 0) + 1
        dupes = [c for c, cnt in seen.items() if cnt > 1]
        if dupes:
            return self._fail(
                f"Duplicate column names found: {dupes}",
                affected_items=[str(d) for d in dupes],
                suggestion="Rename duplicate sample columns to unique identifiers.",
            )
        return self._pass("No duplicate column names found.")


class NonNumericRule(BaseRule):
    rule_id = "FMT-005"
    category = "format"
    severity = "ERROR"
    description = "Count matrix values must be numeric (no text, NA, or empty cells)."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix could not be parsed.")
        df = context.count_matrix
        # Try to coerce to numeric
        numeric_df = df.apply(pd.to_numeric, errors="coerce")
        na_mask = numeric_df.isna() & df.notna()
        non_numeric_count = int(na_mask.values.sum())
        null_count = int(numeric_df.isna().values.sum())

        issues = []
        if non_numeric_count > 0:
            issues.append(f"{non_numeric_count} non-numeric value(s) detected")
        if null_count > 0:
            issues.append(f"{null_count} missing/NA value(s) detected")

        if issues:
            return self._fail(
                "; ".join(issues) + " in count matrix.",
                suggestion=(
                    "Replace non-numeric and missing values with 0 or remove affected rows/columns."
                ),
                details={"non_numeric_count": non_numeric_count, "null_count": null_count},
            )
        return self._pass("All count matrix values are numeric.")


class NegativeCountRule(BaseRule):
    rule_id = "FMT-006"
    category = "format"
    severity = "ERROR"
    description = "Raw count values must be non-negative integers."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix could not be parsed.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce")
        neg_mask = df < 0
        neg_count = int(neg_mask.values.sum())
        if neg_count > 0:
            # Collect gene IDs with negative values
            affected = list(
                df.index[neg_mask.any(axis=1)].astype(str)[:10]
            )
            return self._fail(
                f"{neg_count} negative value(s) found. Raw counts must be ≥ 0.",
                affected_items=affected,
                suggestion="Remove or correct rows with negative counts.",
            )
        return self._pass("No negative count values detected.")


class WhitespaceNameRule(BaseRule):
    rule_id = "FMT-007"
    category = "format"
    severity = "WARNING"
    description = "Sample and gene names should not contain leading/trailing whitespace or special characters."

    _SPECIAL = re.compile(r'[^\w\-.]')

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix could not be parsed.")
        df = context.count_matrix
        issues = []
        for col in df.columns:
            sc = str(col)
            if sc != sc.strip():
                issues.append(f"Column '{sc}' has leading/trailing whitespace")
            elif self._SPECIAL.search(sc):
                issues.append(f"Column '{sc}' contains special characters")
        for idx in df.index[:500]:
            si = str(idx)
            if si != si.strip():
                issues.append(f"Gene ID '{si}' has leading/trailing whitespace")

        if issues:
            return self._fail(
                f"{len(issues)} name(s) with whitespace or special characters.",
                affected_items=issues[:10],
                suggestion="Strip whitespace and replace special characters with underscores.",
            )
        return self._pass("No whitespace or special character issues in names.")
