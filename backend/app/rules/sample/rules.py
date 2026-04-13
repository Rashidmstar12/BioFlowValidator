"""Sample-level validation rules."""
from __future__ import annotations

import numpy as np
import pandas as pd

from app.models.context import ValidationContext
from app.models.rule_result import RuleResult
from app.rules.base import BaseRule


class SampleMatchRule(BaseRule):
    rule_id = "SMP-001"
    category = "sample"
    severity = "ERROR"
    description = "Sample IDs in count matrix must exactly match those in metadata."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        if context.metadata is None:
            return self._skip("Metadata not provided; sample match check skipped.")

        matrix_samples = set(str(c) for c in context.count_matrix.columns)
        meta_samples = set(str(i) for i in context.metadata.index)

        only_matrix = sorted(matrix_samples - meta_samples)
        only_meta = sorted(meta_samples - matrix_samples)

        if only_matrix or only_meta:
            affected = (
                [f"In matrix only: {s}" for s in only_matrix[:5]]
                + [f"In metadata only: {s}" for s in only_meta[:5]]
            )
            return self._fail(
                f"{len(only_matrix)} sample(s) in matrix not in metadata; "
                f"{len(only_meta)} sample(s) in metadata not in matrix.",
                affected_items=affected,
                suggestion=(
                    "Ensure sample IDs are identical (case-sensitive) in both files. "
                    "Check for typos, trailing spaces, or missing samples."
                ),
                details={
                    "only_in_matrix": only_matrix,
                    "only_in_metadata": only_meta,
                },
            )
        return self._pass(
            f"All {len(matrix_samples)} sample IDs match between count matrix and metadata."
        )


class SampleOrderRule(BaseRule):
    rule_id = "SMP-002"
    category = "sample"
    severity = "WARNING"
    description = "Sample order in count matrix should match metadata row order."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None or context.metadata is None:
            return self._skip("Count matrix or metadata not available.")

        matrix_order = [str(c) for c in context.count_matrix.columns]
        meta_order = [str(i) for i in context.metadata.index]

        # Only compare common samples
        common = [s for s in matrix_order if s in set(meta_order)]
        meta_common = [s for s in meta_order if s in set(matrix_order)]

        if common != meta_common:
            return self._fail(
                "Sample order differs between count matrix and metadata. "
                "Some tools (e.g. DESeq2) require identical ordering.",
                suggestion=(
                    "Reorder metadata rows to match the column order of the count matrix, "
                    "or reorder count matrix columns to match metadata."
                ),
            )
        return self._pass("Sample order is consistent between count matrix and metadata.")


class DuplicateSampleRule(BaseRule):
    rule_id = "SMP-003"
    category = "sample"
    severity = "ERROR"
    description = "Sample IDs must be unique in both count matrix and metadata."

    def run(self, context: ValidationContext) -> RuleResult:
        issues = []
        if context.count_matrix is not None:
            cols = [str(c) for c in context.count_matrix.columns]
            seen: dict = {}
            for c in cols:
                seen[c] = seen.get(c, 0) + 1
            dupes = [c for c, n in seen.items() if n > 1]
            if dupes:
                issues.append(f"Duplicate columns in matrix: {dupes[:5]}")

        if context.metadata is not None:
            idx = [str(i) for i in context.metadata.index]
            seen2: dict = {}
            for i in idx:
                seen2[i] = seen2.get(i, 0) + 1
            dupes2 = [i for i, n in seen2.items() if n > 1]
            if dupes2:
                issues.append(f"Duplicate rows in metadata: {dupes2[:5]}")

        if issues:
            return self._fail(
                "Duplicate sample IDs detected.",
                affected_items=issues,
                suggestion="Remove or rename duplicate samples so all IDs are unique.",
            )
        return self._pass("No duplicate sample IDs found.")


class ReplicateCountRule(BaseRule):
    rule_id = "SMP-004"
    category = "sample"
    severity = "ERROR"
    description = "Each condition must have at least 2 biological replicates for DE analysis."

    _CONDITION_COLS = ("condition", "group", "treatment", "genotype", "cell_type")

    def _find_condition_col(self, meta: pd.DataFrame) -> str | None:
        for col in meta.columns:
            if col.lower() in self._CONDITION_COLS:
                return col
        # Fall back to first non-numeric column
        for col in meta.columns:
            if meta[col].dtype == object:
                return col
        return None

    def run(self, context: ValidationContext) -> RuleResult:
        if context.metadata is None:
            return self._skip("Metadata not provided.")
        meta = context.metadata
        cond_col = self._find_condition_col(meta)
        if cond_col is None:
            return self._skip("No condition column found in metadata.")

        counts = meta[cond_col].value_counts()
        under_2 = counts[counts < 2].index.tolist()
        under_3 = counts[(counts >= 2) & (counts < 3)].index.tolist()

        if under_2:
            return self._fail(
                f"Condition(s) {under_2} have fewer than 2 replicates. "
                "DE analysis is statistically invalid with only 1 replicate.",
                affected_items=[str(c) for c in under_2],
                suggestion=(
                    "Include at least 3 biological replicates per condition. "
                    "A minimum of 2 is required, but 3+ is strongly recommended."
                ),
                details={"replicate_counts": counts.to_dict()},
            )
        results = []
        if under_3:
            return RuleResult(
                rule_id=self.rule_id,
                category=self.category,
                severity="WARNING",
                status="FAIL",
                message=(
                    f"Condition(s) {under_3} have only 2 replicates. "
                    "3 or more biological replicates are strongly recommended."
                ),
                affected_items=[str(c) for c in under_3],
                suggestion="Add more biological replicates to increase statistical power.",
                details={"replicate_counts": counts.to_dict()},
            )
        return self._pass(
            f"All conditions have ≥ 3 replicates. Counts: {counts.to_dict()}"
        )


class NearIdenticalSampleRule(BaseRule):
    rule_id = "SMP-005"
    category = "sample"
    severity = "WARNING"
    description = "Detect near-identical sample profiles (possible technical replicate confusion)."

    _CORR_THRESHOLD = 0.9999

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        if df.shape[1] < 2:
            return self._skip("Fewer than 2 samples; correlation check skipped.")

        corr = df.corr()
        issues = []
        cols = list(corr.columns)
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                val = corr.iloc[i, j]
                if val >= self._CORR_THRESHOLD:
                    issues.append(f"'{cols[i]}' and '{cols[j]}' (r={val:.6f})")

        if issues:
            return self._fail(
                f"{len(issues)} pair(s) of samples are near-identical (r ≥ {self._CORR_THRESHOLD}).",
                affected_items=issues[:10],
                suggestion=(
                    "Verify that these are truly independent biological replicates, "
                    "not duplicated technical replicates."
                ),
            )
        return self._pass("No near-identical sample pairs detected.")
