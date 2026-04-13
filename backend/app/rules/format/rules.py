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
        # Heuristic: if columns are sequential integers starting from 0 or 1, the header is
        # likely missing (those are pandas default integer column names). Non-sequential
        # numeric columns (e.g. Entrez IDs) are valid identifiers and must NOT trigger this
        # rule — that would be a false positive on any matrix using Entrez gene IDs as
        # row/column headers.
        try:
            int_cols = [int(c) for c in df.columns]
            first = int_cols[0]
            is_sequential = (
                first in (0, 1)
                and int_cols == list(range(first, first + len(int_cols)))
            )
        except (ValueError, TypeError):
            is_sequential = False

        if is_sequential:
            return self._fail(
                "Column names are sequential integers (0, 1, 2, … or 1, 2, 3, …), "
                "which are pandas default column indices and indicate the header row is missing.",
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


class MatrixOrientationRule(BaseRule):
    """FMT-008 — Detect likely transposed count matrices (samples × genes instead of genes × samples).

    Biological importance: A transposed matrix causes every downstream rule to analyse
    sample names as gene IDs and gene IDs as sample names.  Every rule that fires
    downstream will be a false positive.  This is a catastrophic failure mode that
    makes the tool completely unreliable without this check.

    Heuristic: A typical bulk RNA-seq count matrix has 15,000–25,000 gene rows and
    4–100 sample columns.  If n_columns > n_rows AND n_rows < 500, the matrix is
    almost certainly transposed.  Even targeted gene panels have ≥ 200 genes.
    """

    rule_id = "FMT-008"
    category = "format"
    severity = "WARNING"
    description = (
        "Count matrix orientation should be genes × samples (rows = genes, columns = samples). "
        "If n_columns > n_rows and n_rows < 500, the matrix appears to be transposed."
    )

    _MAX_GENE_ROWS_FOR_TRANSPOSED = 500

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix could not be parsed.")
        df = context.count_matrix
        n_rows, n_cols = df.shape
        if n_cols > n_rows and n_rows < self._MAX_GENE_ROWS_FOR_TRANSPOSED:
            return self._fail(
                f"Count matrix has {n_rows} rows and {n_cols} columns. "
                "Expected orientation is genes × samples (rows = genes). "
                "The matrix appears to be transposed (samples × genes).",
                affected_items=[
                    f"n_rows={n_rows} (expected ~15,000–25,000 for bulk RNA-seq)",
                    f"n_cols={n_cols} (expected ~4–100)",
                ],
                suggestion=(
                    "Transpose the matrix so that rows correspond to genes and columns "
                    "correspond to samples. In Python: df.T"
                ),
                details={"n_rows": n_rows, "n_cols": n_cols},
            )
        return self._pass(
            f"Matrix orientation appears correct: {n_rows} rows (genes) × {n_cols} columns (samples)."
        )
