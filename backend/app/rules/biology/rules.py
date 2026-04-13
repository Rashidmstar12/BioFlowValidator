"""Biological sanity check rules."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from app.models.context import ValidationContext
from app.models.rule_result import RuleResult
from app.rules.base import BaseRule

# Common mitochondrial gene prefix patterns
_MT_PATTERN = re.compile(r'^MT-', re.IGNORECASE)

# Common housekeeping genes (HGNC symbols)
_HOUSEKEEPING_GENES = {
    "ACTB", "GAPDH", "B2M", "HPRT1", "HMBS", "SDHA",
    "TBP", "RPLP0", "YWHAZ", "PPIA",
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
    description = "High mitochondrial gene fraction per sample suggests poor sample quality (relevant for single-cell and bulk RNA-seq)."

    _MT_THRESHOLD = 0.30  # 30%

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        df = context.count_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        mt_mask = pd.Series(
            [bool(_MT_PATTERN.match(str(g))) for g in df.index], index=df.index
        )
        if mt_mask.sum() == 0:
            return self._pass(
                "No mitochondrial genes detected in count matrix "
                "(MT- prefix pattern not found; check may not be applicable)."
            )

        lib_sizes = df.sum(axis=0)
        mt_counts = df[mt_mask].sum(axis=0)
        mt_frac = mt_counts / lib_sizes.replace(0, np.nan)

        high_mt = mt_frac[mt_frac > self._MT_THRESHOLD].dropna()
        if len(high_mt) > 0:
            return self._fail(
                f"{len(high_mt)} sample(s) have mitochondrial fraction > {self._MT_THRESHOLD*100:.0f}%.",
                affected_items=[
                    f"'{s}': {v*100:.1f}% MT" for s, v in high_mt.items()
                ][:10],
                suggestion=(
                    "High MT fraction indicates poor sample quality (cell death, "
                    "cytoplasmic contamination). Consider removing these samples."
                ),
            )
        return self._pass("Mitochondrial gene fractions are within acceptable range.")


class HousekeepingGeneRule(BaseRule):
    rule_id = "BIO-006"
    category = "biology"
    severity = "WARNING"
    description = "Common housekeeping genes should be present if gene symbols are used."

    _ENSEMBL_PATTERN = re.compile(r'^ENSG\d+', re.IGNORECASE)

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        gene_ids = {str(g).upper() for g in context.count_matrix.index}
        # Skip if IDs look like Ensembl IDs
        ensembl_like = sum(1 for g in gene_ids if self._ENSEMBL_PATTERN.match(g))
        if ensembl_like > len(gene_ids) * 0.5:
            return self._skip(
                "Gene IDs appear to be Ensembl IDs; housekeeping gene symbol check skipped."
            )
        # Only run if data looks like gene symbols
        symbol_like = sum(
            1 for g in gene_ids
            if re.match(r'^[A-Z][A-Z0-9\-]{1,}$', g)
        )
        if symbol_like < len(gene_ids) * 0.5:
            return self._skip(
                "Gene IDs do not appear to be gene symbols; housekeeping gene check skipped."
            )

        missing = _HOUSEKEEPING_GENES - gene_ids
        if len(missing) >= len(_HOUSEKEEPING_GENES) * 0.8:
            return self._fail(
                f"{len(missing)} of {len(_HOUSEKEEPING_GENES)} common housekeeping genes absent: {sorted(missing)}.",
                affected_items=sorted(missing),
                suggestion=(
                    "Absence of most housekeeping genes is suspicious. Verify gene ID namespace "
                    "and that the count matrix includes expected genes."
                ),
            )
        return self._pass(
            f"Most housekeeping genes present. Missing: {sorted(missing) or 'none'}."
        )
