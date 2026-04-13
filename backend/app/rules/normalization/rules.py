"""Normalization / statistical validation rules."""
from __future__ import annotations

import numpy as np
import pandas as pd

from app.models.context import ValidationContext
from app.models.rule_result import RuleResult
from app.rules.base import BaseRule


class NonIntegerCountRule(BaseRule):
    rule_id = "NRM-001"
    category = "normalization"
    severity = "ERROR"
    description = "Count matrix values should be non-negative integers (raw counts). Non-integer values suggest pre-normalization."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce")
        # Drop NA
        flat = df.values.flatten()
        flat = flat[~np.isnan(flat)]
        if len(flat) == 0:
            return self._skip("No numeric values to check.")

        non_int = np.sum(flat != np.floor(flat))
        frac = non_int / len(flat)

        if frac > 0.01:  # > 1% non-integer → likely normalised
            return self._fail(
                f"{non_int} ({frac*100:.1f}%) values are non-integer. "
                "This suggests the matrix contains normalised values (RPKM/TPM/CPM) "
                "rather than raw integer counts.",
                suggestion=(
                    "Supply raw integer counts as produced by featureCounts, HTSeq, "
                    "STAR, or Salmon/tximport. Normalisation is applied internally by "
                    "DESeq2/edgeR."
                ),
                details={"non_integer_count": int(non_int), "fraction": float(frac)},
            )
        return self._pass("Count values are integers (consistent with raw counts).")


class LibrarySizeRule(BaseRule):
    rule_id = "NRM-002"
    category = "normalization"
    severity = "WARNING"
    description = "Extreme library size variation (max/min > 10×) between samples may indicate quality issues."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        lib_sizes = df.sum(axis=0)
        if lib_sizes.min() <= 0:
            # Near-zero library will be caught by ZeroLibraryRule
            positive = lib_sizes[lib_sizes > 0]
            if len(positive) < 2:
                return self._skip("Not enough samples with positive library sizes.")
            ratio = positive.max() / positive.min()
        else:
            ratio = lib_sizes.max() / lib_sizes.min()

        if ratio > 10:
            worst = lib_sizes.idxmin()
            return self._fail(
                f"Library size ratio (max/min) = {ratio:.1f}×, exceeding the 10× threshold.",
                affected_items=[f"Smallest library: '{worst}' ({int(lib_sizes[worst]):,} counts)"],
                suggestion=(
                    "Investigate low-count samples for sequencing failures. "
                    "Consider removing samples with library sizes < 500,000 reads."
                ),
                details={
                    "max_library": int(lib_sizes.max()),
                    "min_library": int(lib_sizes.min()),
                    "ratio": float(ratio),
                },
            )
        return self._pass(f"Library size variation is acceptable (ratio = {ratio:.1f}×).")


class ZeroLibraryRule(BaseRule):
    rule_id = "NRM-003"
    category = "normalization"
    severity = "ERROR"
    description = "Samples with near-zero total counts indicate failed libraries."

    _THRESHOLD = 1000

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        lib_sizes = df.sum(axis=0)
        failed = lib_sizes[lib_sizes < self._THRESHOLD]
        if len(failed) > 0:
            return self._fail(
                f"{len(failed)} sample(s) have total counts < {self._THRESHOLD:,}, suggesting failed libraries.",
                affected_items=[f"'{s}': {int(v)} counts" for s, v in failed.items()],
                suggestion="Remove failed samples before running DE analysis.",
                details={"failed_samples": {str(k): int(v) for k, v in failed.items()}},
            )
        return self._pass(f"All samples have ≥ {self._THRESHOLD:,} total counts.")


class AllZeroGeneRule(BaseRule):
    rule_id = "NRM-004"
    category = "normalization"
    severity = "WARNING"
    description = "Genes with all-zero counts across all samples should be filtered out before DE analysis."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        all_zero_mask = (df == 0).all(axis=1)
        n_zero = int(all_zero_mask.sum())
        n_total = len(df)
        if n_zero > 0:
            pct = n_zero / n_total * 100
            return self._fail(
                f"{n_zero} ({pct:.1f}%) genes have all-zero counts and should be filtered.",
                suggestion=(
                    "Apply low-count filtering before DE analysis. "
                    "A common threshold: keep genes with CPM > 1 in at least n samples "
                    "(where n = size of smallest group)."
                ),
                details={"all_zero_count": n_zero, "total_genes": n_total},
            )
        return self._pass("No all-zero gene rows detected.")


class LowCountDominanceRule(BaseRule):
    rule_id = "NRM-005"
    category = "normalization"
    severity = "WARNING"
    description = "More than 50% of genes with median count < 1 suggests very sparse or low-depth data."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        medians = df.median(axis=1)
        low_count = (medians < 1).sum()
        total = len(medians)
        frac = low_count / total if total > 0 else 0
        if frac > 0.5:
            return self._fail(
                f"{low_count} ({frac*100:.1f}%) genes have median count < 1. "
                "The data may be very sparse or from an unusually low-depth library.",
                suggestion=(
                    "Apply low-count filtering. Consider whether sequencing depth is "
                    "sufficient for DE analysis."
                ),
                details={"low_count_genes": int(low_count), "fraction": float(frac)},
            )
        return self._pass(f"Low-count gene fraction is acceptable ({frac*100:.1f}%).")


class DuplicateRowRule(BaseRule):
    rule_id = "NRM-006"
    category = "normalization"
    severity = "WARNING"
    description = "Detect suspiciously identical gene count rows (possible data duplication)."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        # Only check rows with total > 0 (all-zero duplicates caught by AllZeroGeneRule)
        nonzero = df[df.sum(axis=1) > 0]
        duped = nonzero.duplicated(keep=False)
        n_duped = int(duped.sum())
        if n_duped > 0:
            affected = list(nonzero[duped].index[:10].astype(str))
            return self._fail(
                f"{n_duped} gene rows have identical count profiles, suggesting data duplication.",
                affected_items=affected,
                suggestion="Verify that gene rows are not duplicated; aggregate if necessary.",
            )
        return self._pass("No duplicate gene count rows detected.")
