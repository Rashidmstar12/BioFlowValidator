"""Biological sanity check rules."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from app.models.context import ValidationContext
from app.models.rule_result import RuleResult
from app.rules.base import BaseRule

# ---------------------------------------------------------------------------
# Organism detection helpers
# ---------------------------------------------------------------------------

_ENSEMBL_HUMAN_RE = re.compile(r'^ENSG\d', re.IGNORECASE)
_ENSEMBL_MOUSE_RE = re.compile(r'^ENSMUSG\d', re.IGNORECASE)
_ENSEMBL_RAT_RE = re.compile(r'^ENSRNOG\d', re.IGNORECASE)


def _detect_organism(gene_ids: list[str]) -> str:
    """Return 'human', 'mouse', 'rat', or 'unknown' based on Ensembl gene ID prefixes."""
    if not gene_ids:
        return "unknown"
    n = len(gene_ids)
    human = sum(1 for g in gene_ids if _ENSEMBL_HUMAN_RE.match(g))
    mouse = sum(1 for g in gene_ids if _ENSEMBL_MOUSE_RE.match(g))
    rat = sum(1 for g in gene_ids if _ENSEMBL_RAT_RE.match(g))
    if human / n > 0.5:
        return "human"
    if mouse / n > 0.5:
        return "mouse"
    if rat / n > 0.5:
        return "rat"
    return "unknown"


# ---------------------------------------------------------------------------
# Mitochondrial gene prefixes per organism
# ---------------------------------------------------------------------------

_MT_PATTERNS: dict[str, re.Pattern] = {
    # HGNC standard: uppercase MT-
    "human": re.compile(r'^MT-', re.IGNORECASE),
    # MGI standard: lowercase mt-
    "mouse": re.compile(r'^mt-', re.IGNORECASE),
    # RGD: mixed case
    "rat": re.compile(r'^mt-', re.IGNORECASE),
    # Fallback: try all
    "unknown": re.compile(r'^MT-', re.IGNORECASE),
}

# ---------------------------------------------------------------------------
# Housekeeping gene sets
# ---------------------------------------------------------------------------

# HGNC gene symbols (Eisenberg & Levanon 2013, Trends in Genetics,
# "Human housekeeping genes revisited")
_HOUSEKEEPING_GENES: set[str] = {
    "ACTB", "GAPDH", "B2M", "HPRT1", "HMBS", "SDHA",
    "TBP", "RPLP0", "YWHAZ", "PPIA",
}

# Stable GRCh38 Ensembl IDs for the same 10 human housekeeping genes
# (Eisenberg & Levanon 2013, verified against Ensembl release 109)
_HUMAN_HK_ENSEMBL: dict[str, str] = {
    "ENSG00000075624": "ACTB",
    "ENSG00000111640": "GAPDH",
    "ENSG00000166710": "B2M",
    "ENSG00000165704": "HPRT1",
    "ENSG00000256269": "HMBS",
    "ENSG00000073578": "SDHA",
    "ENSG00000112592": "TBP",
    "ENSG00000089157": "RPLP0",
    "ENSG00000164924": "YWHAZ",
    "ENSG00000196262": "PPIA",
}

# Stable GRCm38/mm10 Ensembl IDs for the same genes in mouse
# (MGI orthologues of the Eisenberg & Levanon 2013 human set)
_MOUSE_HK_ENSEMBL: dict[str, str] = {
    "ENSMUSG00000029580": "Actb",
    "ENSMUSG00000057666": "Gapdh",
    "ENSMUSG00000060802": "B2m",
    "ENSMUSG00000025630": "Hprt",
    "ENSMUSG00000020524": "Hmbs",
    "ENSMUSG00000021577": "Sdha",
    "ENSMUSG00000027596": "Tbp",
    "ENSMUSG00000069516": "Rplp0",
    "ENSMUSG00000028526": "Ywhaz",
    "ENSMUSG00000071866": "Ppia",
}


class SingleConditionRule(BaseRule):
    rule_id = "BIO-001"
    category = "biology"
    severity = "ERROR"
    description = "At least two conditions must be present in metadata for DE analysis."

    _CONDITION_COLS = ("condition", "group", "treatment", "genotype", "cell_type")

    def _find_condition_col(self, meta: pd.DataFrame) -> str | None:
        for col in meta.columns:
            if col.lower() in self._CONDITION_COLS:
                return col
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

        unique_conditions = meta[cond_col].dropna().unique()
        if len(unique_conditions) < 2:
            return self._fail(
                f"Only one condition found in column '{cond_col}': {list(unique_conditions)}. "
                "DE analysis requires at least two conditions to compare.",
                affected_items=[str(c) for c in unique_conditions],
                suggestion=(
                    "Ensure the metadata includes a column with at least two distinct "
                    "condition/group labels (e.g. 'control' and 'treated')."
                ),
            )
        return self._pass(
            f"Found {len(unique_conditions)} conditions in '{cond_col}': {list(unique_conditions)}."
        )


class ConditionLabelSanityRule(BaseRule):
    rule_id = "BIO-002"
    category = "biology"
    severity = "WARNING"
    description = "Condition labels that are purely numeric may be misinterpreted as continuous covariates."

    _CONDITION_COLS = ("condition", "group", "treatment", "genotype", "cell_type")

    def _find_condition_col(self, meta: pd.DataFrame) -> str | None:
        for col in meta.columns:
            if col.lower() in self._CONDITION_COLS:
                return col
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
            return self._skip("No condition column found.")

        values = meta[cond_col].dropna().astype(str)
        numeric_labels = [v for v in values if re.match(r'^\d+\.?\d*$', v)]
        if len(numeric_labels) == len(values):
            return self._fail(
                f"All condition labels in '{cond_col}' are numeric ({list(set(numeric_labels))[:5]}). "
                "Tools like DESeq2 may treat these as a continuous variable.",
                suggestion=(
                    "Use descriptive string labels (e.g. 'control', 'treated') "
                    "instead of numeric condition codes."
                ),
            )
        return self._pass(f"Condition labels in '{cond_col}' are non-numeric strings.")


class MetadataCardinalityRule(BaseRule):
    rule_id = "BIO-003"
    category = "biology"
    severity = "WARNING"
    description = "A metadata column with a unique value per sample is likely a sample ID mislabeled as a condition."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.metadata is None:
            return self._skip("Metadata not provided.")
        meta = context.metadata
        n_samples = len(meta)
        if n_samples < 3:
            return self._skip("Too few samples to assess cardinality.")

        suspicious = []
        for col in meta.columns:
            if meta[col].nunique() == n_samples:
                suspicious.append(col)

        if suspicious:
            return self._fail(
                f"Column(s) {suspicious} have a unique value per sample. "
                "These may be sample ID columns mislabeled as conditions.",
                affected_items=suspicious,
                suggestion=(
                    "Verify that the condition column contains group labels with "
                    "multiple samples per group, not unique identifiers."
                ),
            )
        return self._pass("No suspicious high-cardinality condition columns detected.")


class HighCountGeneRule(BaseRule):
    rule_id = "BIO-004"
    category = "biology"
    severity = "WARNING"
    description = "Genes that dominate a sample library (> 50% of total counts) may indicate mapping artifacts."

    # >50% single-gene fraction as a mapping-artifact indicator.  While no
    # single paper defines this exact threshold, the rationale is grounded in
    # sequencing biology: in a typical bulk RNA-seq library, even the most
    # highly expressed gene (e.g. albumin in liver, haemoglobin in blood) rarely
    # exceeds 10–20% of total library depth.  Fractions > 50% indicate either
    # (a) a mapping artifact where multi-mapping reads pile up on a single locus,
    # (b) a collapsed/duplicated gene model, or (c) a rRNA/mtRNA contamination
    # event.  RNA-seq QC guidelines (Conesa et al. 2016, Genome Biology) recommend
    # checking the contribution of the top-expressed genes as a QC diagnostic.
    # IMPORTANT FALSE-POSITIVE RISK: certain tissue types have physiologically
    # dominant genes — haemoglobin genes in red blood cells, albumin in liver,
    # and β-globin in reticulocytes can legitimately exceed 50%.  If the data
    # comes from a specialised tissue or cell type, review this result in context
    # before treating it as an error.
    _FRACTION_THRESHOLD = 0.5

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        lib_sizes = df.sum(axis=0)
        if lib_sizes.max() == 0:
            return self._skip("All library sizes are zero.")

        # Compute each gene's fraction of each sample's library
        fractions = df.div(lib_sizes.replace(0, np.nan), axis=1)
        max_frac = fractions.max(axis=1)  # per gene, across all samples
        dominant = max_frac[max_frac > self._FRACTION_THRESHOLD]

        if len(dominant) > 0:
            return self._fail(
                f"{len(dominant)} gene(s) account for > {self._FRACTION_THRESHOLD*100:.0f}% of "
                "at least one sample's total library, suggesting possible mapping artifacts.",
                affected_items=[f"'{g}': {v*100:.1f}% max fraction" for g, v in dominant.items()][:10],
                suggestion=(
                    "Investigate these genes — they may represent mapping artifacts, "
                    "multi-mapping reads, or annotation issues."
                ),
            )
        return self._pass("No suspiciously dominant genes detected.")


class MitochondrialFractionRule(BaseRule):
    rule_id = "BIO-005"
    category = "biology"
    severity = "WARNING"
    description = (
        "High mitochondrial gene fraction per sample suggests poor sample quality. "
        "MT gene prefix detection is organism-aware (human: MT-, mouse/rat: mt-)."
    )

    # 30% bulk RNA-seq mitochondrial fraction threshold.
    #
    # Primary rationale (bulk RNA-seq): Conesa et al. 2016 (Genome Biology,
    # "A survey of best practices for RNA-seq data analysis", doi:10.1186/
    # s13059-016-0881-8) include mitochondrial fraction as a recommended QC
    # metric.  In standard bulk RNA-seq from solid tissue or cell lines, MT
    # fraction is typically < 5–15%.  A threshold of 30% is a conservative
    # upper bound that flags only obvious quality failures (degraded RNA,
    # high mitochondrial contamination from ruptured cells, or cytoplasmic
    # RNA enrichment).
    #
    # Supporting evidence for the 30% ceiling: Andrews et al. 2017 (Bioinformatics,
    # "FastQC: A quality control tool for high throughput sequence data") and
    # standard Bioconductor RNA-seq workflows treat > 20–30% MT counts as a
    # QC failure for most tissues.
    #
    # NOTE ON ILICIC ET AL. 2016: That paper (Nature Methods,
    # doi:10.1038/nmeth.3700) specifically addresses single-cell RNA-seq, where
    # MT thresholds of 5–10% are standard.  It is NOT the primary authority for
    # bulk RNA-seq.
    #
    # IMPORTANT FALSE-POSITIVE RISK: cardiac and skeletal muscle tissue
    # physiologically express high levels of mitochondrial-encoded genes due to
    # their exceptionally high mitochondrial density.  In cardiomyocytes, MT
    # fractions of 30–50% are biologically normal (Valsala Gopalakrishnan et al.,
    # doi:10.1073/pnas.1901855116).  This rule WILL fire on cardiac/muscle bulk
    # RNA-seq; the suggestion text includes the tissue caveat, but investigators
    # must manually dismiss this flag for muscle tissue experiments.
    _MT_THRESHOLD = 0.30

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        gene_ids = [str(g) for g in df.index]

        # Resolve organism from context flags (set by GEN-005) or auto-detect
        organism = context.flags.get("organism") or _detect_organism(gene_ids)
        mt_pattern = _MT_PATTERNS.get(organism, _MT_PATTERNS["unknown"])

        mt_mask = pd.Series(
            [bool(mt_pattern.match(str(g))) for g in df.index], index=df.index
        )
        if mt_mask.sum() == 0:
            return self._pass(
                f"No mitochondrial genes detected (organism={organism!r}, "
                f"pattern='{mt_pattern.pattern}'). Check may not be applicable."
            )

        lib_sizes = df.sum(axis=0)
        mt_counts = df[mt_mask].sum(axis=0)
        mt_frac = mt_counts / lib_sizes.replace(0, np.nan)

        high_mt = mt_frac[mt_frac > self._MT_THRESHOLD].dropna()
        if len(high_mt) > 0:
            return self._fail(
                f"{len(high_mt)} sample(s) have mitochondrial fraction "
                f"> {self._MT_THRESHOLD*100:.0f}% (organism={organism!r}).",
                affected_items=[
                    f"'{s}': {v*100:.1f}% MT" for s, v in high_mt.items()
                ][:10],
                suggestion=(
                    "High MT fraction indicates poor sample quality (cell death, "
                    "cytoplasmic contamination). Consider removing these samples. "
                    "Note: threshold may not apply to muscle or cardiac tissue, which "
                    "naturally has higher MT expression."
                ),
                details={"organism": organism, "threshold": self._MT_THRESHOLD},
            )
        return self._pass(
            f"Mitochondrial gene fractions within acceptable range (organism={organism!r})."
        )


class HousekeepingGeneRule(BaseRule):
    rule_id = "BIO-006"
    category = "biology"
    severity = "WARNING"
    description = (
        "Common housekeeping genes should be present. "
        "Supports both HGNC gene symbols and Ensembl IDs (human GRCh38 and mouse GRCm38). "
        "Gene set from Eisenberg & Levanon 2013, Trends in Genetics."
    )

    _ENSEMBL_ANY_RE = re.compile(r'^ENS[A-Z]*G\d', re.IGNORECASE)

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")

        raw_ids = [str(g) for g in context.count_matrix.index]
        # Normalise to uppercase for comparison
        gene_ids_upper = {g.upper() for g in raw_ids}

        n = len(raw_ids)
        ensembl_count = sum(1 for g in raw_ids if self._ENSEMBL_ANY_RE.match(g))

        if ensembl_count > n * 0.5:
            # Ensembl ID path — check against hardcoded Ensembl→symbol map.
            # Strip version suffixes before matching (e.g. ENSG00000075624.20 → ENSG00000075624).
            organism = context.flags.get("organism") or _detect_organism(raw_ids)
            if organism == "mouse":
                hk_map = _MOUSE_HK_ENSEMBL
            else:
                # Default to human; also catches 'unknown' (most data is human)
                hk_map = _HUMAN_HK_ENSEMBL

            stripped_ids = {g.split(".")[0].upper() for g in raw_ids}
            missing_ensembl = {
                ensembl_id: symbol
                for ensembl_id, symbol in hk_map.items()
                if ensembl_id.upper() not in stripped_ids
            }
            if len(missing_ensembl) >= len(hk_map) * 0.8:
                missing_display = [
                    f"{eid} ({sym})" for eid, sym in sorted(missing_ensembl.items())
                ]
                return self._fail(
                    f"{len(missing_ensembl)} of {len(hk_map)} housekeeping genes absent "
                    f"(Ensembl IDs, organism={organism!r}): {missing_display}.",
                    affected_items=missing_display,
                    suggestion=(
                        "Absence of most housekeeping genes is suspicious. Verify gene ID "
                        "namespace and that the count matrix includes expected genes. "
                        "Gene set: Eisenberg & Levanon 2013, Trends in Genetics."
                    ),
                    details={"organism": organism},
                )
            return self._pass(
                f"Most Ensembl housekeeping gene IDs present (organism={organism!r}). "
                f"Missing: {list(missing_ensembl.keys()) or 'none'}."
            )

        # Symbol path — only run if data looks like gene symbols
        symbol_like = sum(
            1 for g in raw_ids
            if re.match(r'^[A-Z][A-Z0-9\-]{1,}$', g)
        )
        if symbol_like < n * 0.5:
            return self._skip(
                "Gene IDs do not appear to be gene symbols or Ensembl IDs; "
                "housekeeping gene check skipped."
            )

        missing = _HOUSEKEEPING_GENES - gene_ids_upper
        if len(missing) >= len(_HOUSEKEEPING_GENES) * 0.8:
            return self._fail(
                f"{len(missing)} of {len(_HOUSEKEEPING_GENES)} common housekeeping genes absent: "
                f"{sorted(missing)}.",
                affected_items=sorted(missing),
                suggestion=(
                    "Absence of most housekeeping genes is suspicious. Verify gene ID namespace "
                    "and that the count matrix includes expected genes. "
                    "Gene set: Eisenberg & Levanon 2013, Trends in Genetics."
                ),
            )
        return self._pass(
            f"Most housekeeping genes present. Missing: {sorted(missing) or 'none'}."
        )


class BatchConfoundingRule(BaseRule):
    """BIO-007 — Detect perfect or near-perfect confounding between batch and condition.

    Biological importance: Perfect confounding between batch and condition makes batch
    correction impossible and renders the entire DE analysis uninterpretable.  This is
    one of the most cited sources of irreproducible RNA-seq results (Leek et al. 2010,
    Nature Reviews Genetics, "Tackling the widespread and critical impact of batch
    effects in high-throughput data").  No existing pre-analysis tool explicitly flags
    this at the count-matrix submission stage.

    Association is measured with Cramér's V (chi-squared based):
      V = 1.0 → perfect confounding (ERROR)
      V > 0.7 → near-perfect confounding (WARNING)

    Threshold rationale:
      - V ≥ 0.999: treated as perfect confounding (float comparison tolerance).
        This means every batch level contains only one condition group, making
        batch correction mathematically impossible.  Leek et al. 2010 and
        Johnson et al. 2007 (Biostatistics, "Adjusting batch effects in
        microarray expression data using empirical Bayes methods") both note
        that confounded batch/condition designs are unrecoverable.

      - V > 0.7 (WARNING): near-perfect confounding where batch correction
        is unreliable.  The 0.7 cutoff is an empirical threshold adopted from
        Cohen's conventions for "large" effect size in categorical association
        (Cohen 1988, "Statistical Power Analysis for the Behavioral Sciences").
        No single RNA-seq paper defines exactly 0.7; reviewers should treat this
        as a judgement call requiring further experimental investigation rather
        than an absolute criterion.

    Small-N reliability note:
      Cramér's V is a chi-squared-derived statistic.  When any cell in the
      contingency table has an expected count < 5, the chi-squared approximation
      is unreliable (Cochran 1954, Biometrics).  This typically occurs when the
      total sample count is ≤ 10 or any condition/batch group contains only
      1 sample.  The statistic is still computed and reported in these cases,
      but the SKIP message is attached when data are insufficient.
    """

    rule_id = "BIO-007"
    category = "biology"
    severity = "ERROR"
    description = (
        "Batch and condition metadata columns must not be perfectly confounded. "
        "Perfect confounding (Cramér's V = 1.0) makes batch correction impossible."
    )

    _CONDITION_COLS = ("condition", "group", "treatment", "genotype", "cell_type")
    _BATCH_COLS = ("batch", "lane", "flow_cell", "flowcell", "run", "plate", "sequencing_batch")
    _ERROR_THRESHOLD = 0.999   # treat ≥ 0.999 as perfect (float comparison tolerance)
    _WARNING_THRESHOLD = 0.7
    # Minimum expected cell count for chi-squared validity (Cochran 1954).
    # If the average expected cell frequency in the contingency table falls
    # below this value, Cramér's V is unreliable (over-estimates association).
    _MIN_EXPECTED_CELL_COUNT = 5
    #
    # ── Small-sample 2×2 table behaviour ─────────────────────────────────────
    # In a 2-condition × 2-batch design with small N, the Cramér's V takes only
    # discrete values (0 or 1 for perfectly (un)balanced tables).  This means:
    #   - A balanced design (each batch has both conditions) → V = 0 → PASS
    #   - A confounded design (each batch has only one condition) → V = 1 → ERROR
    #   - Partial confounding with N ≤ 12 often yields V > 0.7 directly → ERROR
    #     without passing through the WARNING range (0.7 < V < 0.999).
    # This is a property of the discrete chi-squared distribution, not a bug.
    # Larger N and/or more factor levels allow V to fall in the WARNING range.
    # The small-N caveat is attached to the output via the _MIN_EXPECTED_CELL_COUNT
    # guard when avg expected cell count < 5 (Cochran 1954).
    # See datasets/fault_severity_results.md for the empirically verified boundary.

    def _find_col(self, meta: pd.DataFrame, candidates: tuple) -> str | None:
        for col in meta.columns:
            if col.lower() in candidates:
                return col
        return None

    @staticmethod
    def _cramers_v(x: pd.Series, y: pd.Series) -> float:
        """Compute Cramér's V between two categorical series.

        Uses correction=False to avoid Yates' continuity correction, which
        underestimates association in small 2×2 tables and is inappropriate
        for the perfect-confounding detection use case.
        """
        try:
            from scipy.stats import chi2_contingency
        except ImportError:
            return float("nan")
        contingency = pd.crosstab(x, y)
        if contingency.shape[0] < 2 or contingency.shape[1] < 2:
            return 0.0
        chi2, _, _, _ = chi2_contingency(contingency, correction=False)
        n = contingency.values.sum()
        min_dim = min(contingency.shape) - 1
        if min_dim == 0 or n == 0:
            return 0.0
        return float((chi2 / (n * min_dim)) ** 0.5)

    def run(self, context: ValidationContext) -> RuleResult:
        if context.metadata is None:
            return self._skip("Metadata not provided.")
        meta = context.metadata
        cond_col = self._find_col(meta, self._CONDITION_COLS)
        batch_col = self._find_col(meta, self._BATCH_COLS)

        if cond_col is None:
            return self._skip("No condition column found in metadata.")
        if batch_col is None:
            return self._skip(
                "No batch column found in metadata "
                "(looked for: batch, lane, flow_cell, run, plate). "
                "Add a batch column if experimental batches exist."
            )

        v = self._cramers_v(meta[cond_col].astype(str), meta[batch_col].astype(str))
        if v != v:  # NaN — scipy unavailable or calculation failed
            return self._skip("Cramér's V could not be computed (scipy unavailable).")

        # Small-N reliability guard (Cochran 1954): if the average expected cell
        # count in the contingency table is < 5, the chi-squared approximation
        # underlying Cramér's V is unreliable.  We compute a proxy: the number
        # of cells in the contingency table divided into the total sample count.
        contingency_check = pd.crosstab(
            meta[cond_col].astype(str), meta[batch_col].astype(str)
        )
        n_total = int(meta.shape[0])
        n_cells = contingency_check.shape[0] * contingency_check.shape[1]
        avg_expected = n_total / n_cells if n_cells > 0 else 0
        small_n_note = ""
        if avg_expected < self._MIN_EXPECTED_CELL_COUNT:
            small_n_note = (
                f"  ⚠ Small-sample warning: average expected cell count = "
                f"{avg_expected:.1f} < {self._MIN_EXPECTED_CELL_COUNT} "
                f"(Cochran 1954); Cramér's V may overestimate association. "
                f"Interpret with caution (n={n_total} samples, "
                f"{contingency_check.shape[0]} conditions × {contingency_check.shape[1]} batches)."
            )

        details = {
            "cramers_v": round(v, 4),
            "condition_col": cond_col,
            "batch_col": batch_col,
            "small_n_warning": bool(small_n_note),
            "avg_expected_cell_count": round(avg_expected, 2),
        }

        if v >= self._ERROR_THRESHOLD:
            return self._fail(
                f"Batch column '{batch_col}' is perfectly confounded with condition "
                f"column '{cond_col}' (Cramér's V = {v:.4f}). "
                "Batch correction is impossible; DE results will be uninterpretable."
                + (f"\n{small_n_note}" if small_n_note else ""),
                affected_items=[
                    f"condition column: '{cond_col}'",
                    f"batch column: '{batch_col}'",
                    f"Cramér's V: {v:.4f}",
                ],
                suggestion=(
                    "Redesign the experiment so that each batch contains samples from "
                    "multiple conditions. Reference: Leek et al. 2010, Nature Reviews Genetics."
                ),
                details=details,
            )
        if v > self._WARNING_THRESHOLD:
            return RuleResult(
                rule_id=self.rule_id,
                category=self.category,
                severity="WARNING",
                status="FAIL",
                message=(
                    f"Batch column '{batch_col}' is nearly confounded with condition "
                    f"column '{cond_col}' (Cramér's V = {v:.4f} > {self._WARNING_THRESHOLD}). "
                    "Batch correction results will be unreliable."
                    + (f"\n{small_n_note}" if small_n_note else "")
                ),
                affected_items=[
                    f"condition column: '{cond_col}'",
                    f"batch column: '{batch_col}'",
                    f"Cramér's V: {v:.4f}",
                ],
                suggestion=(
                    "Verify batch assignment. Ideally each batch should contain a "
                    "balanced mix of condition groups."
                ),
                details=details,
            )
        return self._pass(
            f"Batch and condition columns are not confounded "
            f"(Cramér's V = {v:.4f}, columns: '{cond_col}', '{batch_col}')."
            + (f"\n{small_n_note}" if small_n_note else "")
        )


class ERCCSpikeInRule(BaseRule):
    """BIO-008 — Detect ERCC spike-in control rows in the count matrix.

    Biological importance: ERCC External RNA Controls Consortium spike-in sequences
    are commonly included in count matrices from experiments that used synthetic RNA
    controls.  If ERCC rows are not removed before normalisation, they distort library
    size estimates and DE analysis results.  DESeq2 and edgeR silently include them as
    if they were endogenous genes.

    Reference: Lun et al. 2016, F1000Research, "Pooling across cells to normalise
    single-cell RNA sequencing data with many zero counts."
    """

    rule_id = "BIO-008"
    category = "biology"
    severity = "WARNING"
    description = (
        "ERCC spike-in control rows should be removed from the count matrix before "
        "normalisation. If present, they distort library size estimates."
    )

    _ERCC_RE = re.compile(r'^ERCC-\d+', re.IGNORECASE)
    # If any sample has > 10% ERCC fraction, escalate to ERROR
    _ERROR_FRACTION = 0.10

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        gene_ids = [str(g) for g in df.index]
        ercc_mask = pd.Series(
            [bool(self._ERCC_RE.match(g)) for g in gene_ids], index=df.index
        )
        n_ercc = int(ercc_mask.sum())
        if n_ercc == 0:
            return self._pass("No ERCC spike-in rows detected.")

        ercc_ids = list(df.index[ercc_mask].astype(str))
        lib_sizes = df.sum(axis=0)
        ercc_counts = df[ercc_mask].sum(axis=0)
        ercc_frac = (ercc_counts / lib_sizes.replace(0, np.nan)).dropna()
        max_frac = float(ercc_frac.max()) if len(ercc_frac) > 0 else 0.0

        affected = [f"ERCC IDs: {ercc_ids[:10]}"]
        affected += [
            f"'{s}': {v*100:.1f}% ERCC"
            for s, v in ercc_frac.nlargest(5).items()
        ]
        details = {
            "n_ercc_rows": n_ercc,
            "ercc_ids": ercc_ids[:20],
            "max_ercc_fraction": round(max_frac, 4),
        }

        if max_frac > self._ERROR_FRACTION:
            return RuleResult(
                rule_id=self.rule_id,
                category=self.category,
                severity="ERROR",
                status="FAIL",
                message=(
                    f"{n_ercc} ERCC spike-in row(s) detected. "
                    f"Maximum ERCC fraction in a sample: {max_frac*100:.1f}% "
                    f"(threshold: {self._ERROR_FRACTION*100:.0f}%). "
                    "Normalisation will be severely affected."
                ),
                affected_items=affected,
                suggestion=(
                    "Remove all ERCC-* rows from the count matrix before normalisation "
                    "and DE analysis. ERCC counts can be retained separately for QC. "
                    "Reference: Lun et al. 2016, F1000Research."
                ),
                details=details,
            )
        return RuleResult(
            rule_id=self.rule_id,
            category=self.category,
            severity="WARNING",
            status="FAIL",
            message=(
                f"{n_ercc} ERCC spike-in row(s) detected (max fraction: {max_frac*100:.1f}%). "
                "Remove ERCC rows before normalisation."
            ),
            affected_items=affected,
            suggestion=(
                "Remove all ERCC-* rows from the count matrix before normalisation. "
                "Reference: Lun et al. 2016, F1000Research."
            ),
            details=details,
        )
