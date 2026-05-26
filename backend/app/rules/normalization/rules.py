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

        if frac > 0.01:  # > 1% non-integer → likely normalised or expected counts
            col_sums = df.sum(axis=0)
            # Check if all sums are very close to 1,000,000 (typical of TPM/CPM)
            is_cpm_tpm = all(9.9e5 <= s <= 1.01e6 for s in col_sums)
            # Check if average library size is very small (typical of FPKM or small subsets)
            mean_lib_size = col_sums.mean()
            # If average library size is small (< 500,000)
            is_fpkm = mean_lib_size < 500000

            if is_cpm_tpm or is_fpkm:
                # Pre-normalized data is a fatal ERROR
                return self._fail(
                    f"{non_int} ({frac*100:.1f}%) values are non-integer, and library depth indicates pre-normalized data. "
                    "This suggests the matrix contains normalised values (RPKM/TPM/CPM) "
                    "rather than raw integer counts.",
                    suggestion=(
                        "Supply raw integer counts as produced by featureCounts, HTSeq, "
                        "STAR, or Salmon/tximport. Normalisation is applied internally by "
                        "DESeq2/edgeR."
                    ),
                    details={"non_integer_count": int(non_int), "fraction": float(frac), "mean_library_size": float(mean_lib_size)},
                )
            else:
                # Raw expected counts are only a WARNING
                return RuleResult(
                    rule_id=self.rule_id,
                    category=self.category,
                    severity="WARNING",
                    status="FAIL",
                    message=(
                        f"{non_int} ({frac*100:.1f}%) values are non-integer, but library depth suggests expected counts. "
                        "This is common for RSEM, Kallisto, or Salmon expected counts. "
                        "Note that while some downstream DE tools (e.g. edgeR) accept fractional counts, "
                        "others (e.g. DESeq2) require rounding or import via tximport."
                    ),
                    suggestion=(
                        "If using DESeq2, ensure expected counts are rounded to integers or "
                        "imported using tximport. edgeR and limma-voom can handle fractional expected counts directly."
                    ),
                    details={"non_integer_count": int(non_int), "fraction": float(frac), "mean_library_size": float(mean_lib_size)},
                )
        return self._pass("Count values are integers (consistent with raw counts).")


class LibrarySizeRule(BaseRule):
    rule_id = "NRM-002"
    category = "normalization"
    severity = "WARNING"
    description = "Extreme library size variation (max/min > 10×) between samples may indicate quality issues."

    # 10× max/min library size threshold.  Conesa et al. 2016 (Genome Biology,
    # "A survey of best practices for RNA-seq data analysis") recommend
    # inspecting library size distributions as a primary QC step and note that
    # large imbalances introduce biases that TMM/RLE normalization may not fully
    # correct.  Robinson & Oshlack 2010 (Genome Biology, "A scaling normalization
    # method for differential expression analysis of RNA-seq data") show that
    # highly unequal library sizes inflate false-positive rates in DE analysis.
    # The 10× ratio is a widely-used empirical threshold in community workflows
    # (e.g. Bioconductor RNA-seq vignettes); no single paper defines it as a
    # hard cutoff.  Reviewers should note that 5× imbalance can also bias
    # DESeq2 size factors in small experiments — this threshold is conservative.
    #
    # Threshold comparison is STRICT (> 10, not >= 10):
    #   ratio = 10.000 → PASS (exact equality does not trigger the warning)
    #   ratio = 10.001 → FAIL
    # This is intentional to avoid false positives at the exact boundary, which
    # can arise from rounding when counts are multiplied by a scale factor.
    # See datasets/fault_severity_results.md for the detected boundary.
    _RATIO_THRESHOLD = 10

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

        if ratio > self._RATIO_THRESHOLD:
            worst = lib_sizes.idxmin()
            return self._fail(
                f"Library size ratio (max/min) = {ratio:.1f}×, exceeding the "
                f"{self._RATIO_THRESHOLD}× threshold.",
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

    # 1,000-count threshold for a failed library.  While sequencing best
    # practices (e.g. ENCODE RNA-seq standards, 2011; doi:10.1101/gr.136184.111)
    # require > 10 million mapped reads per sample, a total count below 1,000
    # is diagnostic of a catastrophic library failure (no RNA recovered, failed
    # ligation, adapter contamination covering the entire signal).  This is a
    # conservative lower-bound sentinel rather than a quality-sufficiency
    # criterion.  Samples passing this gate may still have inadequate depth;
    # the ENCODE minimum of 30M reads per sample remains the authoritative
    # sufficiency guideline.  Caveat: targeted panels or small spike-in sets
    # can legitimately have total counts below 1,000 — add those assay types
    # to the expected_skips list in your benchmark YAML.
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

    # All-zero gene filtering is recommended in all major DE analysis workflows.
    # Chen et al. 2016 (F1000Research, "From reads to genes to pathways:
    # differential expression analysis of RNA-seq experiments using Rsubread and
    # the edgeR quasi-likelihood pipeline") explicitly filter genes with zero
    # counts in all samples as the first step.  Love et al. 2014 (Genome
    # Biology, "Moderated estimation of fold change and dispersion for RNA-seq
    # data with DESeq2") note that genes not expressed in any sample contribute
    # no information and increase the multiple-testing burden.  False-positive
    # risk: targeted panels or reduced-representation libraries may intentionally
    # omit most genes; the fraction (not just presence) of all-zero genes is the
    # key signal.  A library with 80%+ all-zero genes warrants further review.

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

    # >50% low-count gene fraction as a sparsity indicator.  Robinson & Oshlack
    # 2010 (Genome Biology) and Chen et al. 2016 (F1000Research) recommend
    # retaining only genes with CPM > 1 in at least N samples (where N equals
    # the smallest group size) before DE analysis.  In a typical ~30 M read
    # library this corresponds to a raw count threshold of ~6 per gene per
    # sample.  If the median count is < 1, the gene is rarely detected even at
    # this sequencing depth, and most statistical models will produce unreliable
    # dispersion estimates.  The 50% fraction threshold is an empirical alert
    # that the entire library may be under-sequenced or the wrong file supplied.
    # Caveat: single-cell pseudo-bulk or low-input libraries may legitimately
    # have high sparsity; apply this threshold only to standard bulk RNA-seq.

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
