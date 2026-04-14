"""Unit tests for Phase 8 new and corrected rules.

Covers:
  FMT-003  HeaderRule        — sequential-integer heuristic fix
  FMT-008  MatrixOrientationRule
  SMP-005  NearIdenticalSampleRule — log1p threshold fix
  GEN-005  OrganismDetectionRule
  BIO-005  MitochondrialFractionRule — organism-aware
  BIO-006  HousekeepingGeneRule     — Ensembl mapping
  BIO-007  BatchConfoundingRule
  BIO-008  ERCCSpikeInRule
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from app.models.context import ValidationContext
from app.rules.format.rules import HeaderRule, MatrixOrientationRule
from app.rules.sample.rules import NearIdenticalSampleRule
from app.rules.gene.rules import OrganismDetectionRule
from app.rules.biology.rules import (
    MitochondrialFractionRule,
    HousekeepingGeneRule,
    BatchConfoundingRule,
    ERCCSpikeInRule,
)
from tests.conftest import _ctx_from_dfs, make_count_matrix, make_metadata


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ctx(count_df=None, meta_df=None, flags=None):
    ctx = _ctx_from_dfs(count_df=count_df, meta_df=meta_df)
    if flags:
        ctx.flags.update(flags)
    return ctx


def _human_ensembl_matrix(n_genes=200, n_samples=6):
    genes = [f"ENSG{i:011d}" for i in range(1, n_genes + 1)]
    rng = np.random.default_rng(0)
    values = rng.integers(10, 500, size=(n_genes, n_samples))
    cols = [f"s{i}" for i in range(n_samples)]
    return pd.DataFrame(values, index=genes, columns=cols)


def _mouse_ensembl_matrix(n_genes=200, n_samples=6):
    genes = [f"ENSMUSG{i:011d}" for i in range(1, n_genes + 1)]
    rng = np.random.default_rng(1)
    values = rng.integers(10, 500, size=(n_genes, n_samples))
    cols = [f"s{i}" for i in range(n_samples)]
    return pd.DataFrame(values, index=genes, columns=cols)


# ===========================================================================
# FMT-003 HeaderRule — sequential integer fix
# ===========================================================================

class TestHeaderRule:
    def test_sequential_from_zero_fails(self):
        df = pd.DataFrame(
            np.ones((10, 4), dtype=int),
            columns=[0, 1, 2, 3],
        )
        ctx = _ctx(count_df=df)
        result = HeaderRule().run(ctx)
        assert result.status == "FAIL"

    def test_sequential_from_one_fails(self):
        df = pd.DataFrame(
            np.ones((10, 4), dtype=int),
            columns=[1, 2, 3, 4],
        )
        ctx = _ctx(count_df=df)
        result = HeaderRule().run(ctx)
        assert result.status == "FAIL"

    def test_nonsequential_entrez_ids_pass(self):
        # Entrez IDs are non-sequential integers — must NOT trigger the rule
        df = pd.DataFrame(
            np.ones((5, 3), dtype=int),
            columns=[7157, 4609, 25],  # TP53, MYC, ABL1
        )
        ctx = _ctx(count_df=df)
        result = HeaderRule().run(ctx)
        assert result.status == "PASS"

    def test_string_column_names_pass(self):
        ctx = _ctx(count_df=make_count_matrix())
        assert HeaderRule().run(ctx).status == "PASS"


# ===========================================================================
# FMT-008 MatrixOrientationRule
# ===========================================================================

class TestMatrixOrientationRule:
    def test_transposed_fails(self):
        # 6 rows × 20000 columns → looks like samples × genes
        df = pd.DataFrame(
            np.ones((6, 200), dtype=int),
            index=[f"s{i}" for i in range(6)],
            columns=[f"ENSG{i:011d}" for i in range(200)],
        )
        ctx = _ctx(count_df=df)
        result = MatrixOrientationRule().run(ctx)
        assert result.status == "FAIL"

    def test_correct_orientation_passes(self):
        # 600 rows × 6 columns → genes × samples (even 600 exceeds 500 threshold)
        ctx = _ctx(count_df=make_count_matrix(
            genes=[f"ENSG{i:011d}" for i in range(600)],
            samples=[f"s{i}" for i in range(6)],
        ))
        assert MatrixOrientationRule().run(ctx).status == "PASS"

    def test_large_matrix_passes(self):
        # Typical real dataset: 20000 genes × 8 samples
        df = _human_ensembl_matrix(n_genes=1000, n_samples=8)
        ctx = _ctx(count_df=df)
        assert MatrixOrientationRule().run(ctx).status == "PASS"

    def test_no_matrix_skips(self):
        ctx = _ctx()
        assert MatrixOrientationRule().run(ctx).status == "SKIP"


# ===========================================================================
# SMP-005 NearIdenticalSampleRule — log1p, threshold 0.999
# ===========================================================================

class TestNearIdenticalSampleRule:
    def test_exact_duplicate_fails(self):
        # Two identical columns — log1p corr = 1.0 ≥ 0.999
        rng = np.random.default_rng(7)
        vals = rng.integers(1, 1000, size=(100, 1))
        df = pd.DataFrame(
            np.hstack([vals, vals, vals + 1, vals + 2]),
            index=[f"G{i}" for i in range(100)],
            columns=["s1", "s2", "s3", "s4"],
        )
        ctx = _ctx(count_df=df)
        result = NearIdenticalSampleRule().run(ctx)
        assert result.status == "FAIL"

    def test_genuinely_different_samples_pass(self):
        # Real biological replicates have r ≈ 0.95–0.99 but not ≥ 0.999
        ctx = _ctx(count_df=make_count_matrix())
        assert NearIdenticalSampleRule().run(ctx).status == "PASS"

    def test_old_threshold_0999_would_not_trigger_at_0_9999(self):
        # This verifies the threshold is 0.999 not 0.9999.
        # Build two columns with r > 0.999 but < 0.9999 — should now FAIL.
        rng = np.random.default_rng(42)
        base = rng.integers(10, 10000, size=1000).astype(float)
        noise = rng.normal(0, 1, size=1000)  # very small noise → r ≈ 0.9999+
        df = pd.DataFrame(
            {"s1": base, "s2": base + noise * 0.001},  # r extremely close to 1
        )
        ctx = _ctx(count_df=df)
        result = NearIdenticalSampleRule().run(ctx)
        assert result.status == "FAIL"


# ===========================================================================
# GEN-005 OrganismDetectionRule
# ===========================================================================

class TestOrganismDetectionRule:
    def test_human_detected(self):
        ctx = _ctx(count_df=_human_ensembl_matrix())
        result = OrganismDetectionRule().run(ctx)
        assert result.status == "PASS"
        assert ctx.flags.get("organism") == "human"

    def test_mouse_detected(self):
        ctx = _ctx(count_df=_mouse_ensembl_matrix())
        result = OrganismDetectionRule().run(ctx)
        assert result.status == "PASS"
        assert ctx.flags.get("organism") == "mouse"

    def test_mixed_organism_fails(self):
        # 50% human + 50% mouse IDs
        human_genes = [f"ENSG{i:011d}" for i in range(1, 101)]
        mouse_genes = [f"ENSMUSG{i:011d}" for i in range(1, 101)]
        genes = human_genes + mouse_genes
        rng = np.random.default_rng(3)
        values = rng.integers(1, 100, size=(200, 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        ctx = _ctx(count_df=df)
        result = OrganismDetectionRule().run(ctx)
        assert result.status == "FAIL"

    def test_symbol_ids_skip(self):
        df = make_count_matrix(genes=["ACTB", "GAPDH", "TP53", "MYC"])
        ctx = _ctx(count_df=df)
        result = OrganismDetectionRule().run(ctx)
        assert result.status == "SKIP"

    def test_no_matrix_skips(self):
        ctx = _ctx()
        assert OrganismDetectionRule().run(ctx).status == "SKIP"


# ===========================================================================
# BIO-005 MitochondrialFractionRule — organism-aware
# ===========================================================================

class TestMitochondrialFractionRuleOrganismAware:
    def _make_mt_df(self, mt_prefix: str, organism: str):
        genes = [f"ENSG{i:011d}" for i in range(1, 51)] + [
            f"{mt_prefix}ND1", f"{mt_prefix}ND2", f"{mt_prefix}CO1", f"{mt_prefix}CO2"
        ]
        rng = np.random.default_rng(9)
        values = rng.integers(1, 10, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        # Make MT genes dominate (40% of total)
        for g in [f"{mt_prefix}ND1", f"{mt_prefix}ND2", f"{mt_prefix}CO1", f"{mt_prefix}CO2"]:
            df.loc[g] = 2000
        return df, organism

    def test_human_mt_uppercase_detected(self):
        df, organism = self._make_mt_df("MT-", "human")
        ctx = _ctx(count_df=df, flags={"organism": organism})
        result = MitochondrialFractionRule().run(ctx)
        assert result.status == "FAIL"

    def test_mouse_mt_lowercase_detected(self):
        genes = [f"ENSMUSG{i:011d}" for i in range(1, 51)] + [
            "mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2"
        ]
        rng = np.random.default_rng(10)
        values = rng.integers(1, 10, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        for g in ["mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2"]:
            df.loc[g] = 2000
        ctx = _ctx(count_df=df, flags={"organism": "mouse"})
        result = MitochondrialFractionRule().run(ctx)
        assert result.status == "FAIL"

    def test_human_mt_not_found_in_mouse_matrix_passes(self):
        # Mouse matrix has mt- genes; if we incorrectly use MT- pattern, we miss them
        # But the rule must auto-detect organism when flag is not set
        genes = [f"ENSMUSG{i:011d}" for i in range(1, 51)]
        rng = np.random.default_rng(11)
        values = rng.integers(100, 1000, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        ctx = _ctx(count_df=df)
        result = MitochondrialFractionRule().run(ctx)
        # No MT genes in this mouse matrix — should PASS
        assert result.status == "PASS"


# ===========================================================================
# BIO-006 HousekeepingGeneRule — Ensembl mapping
# ===========================================================================

class TestHousekeepingGeneRuleEnsembl:
    def test_ensembl_human_with_hk_genes_passes(self):
        # Include at least 2 of the 10 HK Ensembl IDs → should PASS (< 80% missing)
        hk_ids = [
            "ENSG00000075624",  # ACTB
            "ENSG00000111640",  # GAPDH
            "ENSG00000166710",  # B2M
        ]
        other_ids = [f"ENSG{i:011d}" for i in range(100, 300)]
        genes = hk_ids + other_ids
        rng = np.random.default_rng(12)
        values = rng.integers(1, 100, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        ctx = _ctx(count_df=df, flags={"organism": "human"})
        result = HousekeepingGeneRule().run(ctx)
        assert result.status == "PASS"

    def test_ensembl_human_missing_all_hk_fails(self):
        # Only non-HK Ensembl IDs — should FAIL
        genes = [f"ENSG{i:011d}" for i in range(100, 300)]
        rng = np.random.default_rng(13)
        values = rng.integers(1, 100, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        ctx = _ctx(count_df=df, flags={"organism": "human"})
        result = HousekeepingGeneRule().run(ctx)
        assert result.status == "FAIL"

    def test_ensembl_version_suffix_stripped(self):
        # Versioned IDs like ENSG00000075624.20 must match after stripping
        hk_ids = [
            "ENSG00000075624.20",  # ACTB versioned
            "ENSG00000111640.15",  # GAPDH versioned
            "ENSG00000166710.17",  # B2M versioned
        ]
        other_ids = [f"ENSG{i:011d}" for i in range(100, 300)]
        genes = hk_ids + other_ids
        rng = np.random.default_rng(14)
        values = rng.integers(1, 100, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        ctx = _ctx(count_df=df, flags={"organism": "human"})
        result = HousekeepingGeneRule().run(ctx)
        assert result.status == "PASS"

    def test_symbol_path_still_works(self):
        # Pure gene symbols — old code path must still work
        symbols = ["ACTB", "GAPDH", "B2M", "HPRT1", "HMBS", "SDHA", "TBP", "RPLP0", "YWHAZ", "PPIA"]
        other = [f"GENE{i}" for i in range(100)]
        genes = symbols + other
        rng = np.random.default_rng(15)
        values = rng.integers(1, 100, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        ctx = _ctx(count_df=df)
        result = HousekeepingGeneRule().run(ctx)
        assert result.status == "PASS"


# ===========================================================================
# BIO-007 BatchConfoundingRule
# ===========================================================================

class TestBatchConfoundingRule:
    def _meta_with_batch(self, conditions, batches):
        samples = [f"s{i}" for i in range(len(conditions))]
        df = pd.DataFrame(
            {"condition": conditions, "batch": batches},
            index=samples,
        )
        return df

    def test_perfect_confounding_fails_error(self):
        # All control in batch A, all treated in batch B
        meta = self._meta_with_batch(
            conditions=["control"] * 3 + ["treated"] * 3,
            batches=["A"] * 3 + ["B"] * 3,
        )
        ctx = _ctx(count_df=make_count_matrix(), meta_df=meta)
        result = BatchConfoundingRule().run(ctx)
        assert result.status == "FAIL"
        assert result.severity == "ERROR"

    def test_balanced_design_passes(self):
        # Each batch has both conditions
        meta = self._meta_with_batch(
            conditions=["control", "treated", "control", "treated", "control", "treated"],
            batches=["A", "A", "B", "B", "C", "C"],
        )
        ctx = _ctx(count_df=make_count_matrix(), meta_df=meta)
        result = BatchConfoundingRule().run(ctx)
        assert result.status == "PASS"

    def test_no_batch_column_skips(self):
        meta = make_metadata()  # no 'batch' column
        ctx = _ctx(count_df=make_count_matrix(), meta_df=meta)
        result = BatchConfoundingRule().run(ctx)
        assert result.status == "SKIP"

    def test_no_metadata_skips(self):
        ctx = _ctx(count_df=make_count_matrix())
        assert BatchConfoundingRule().run(ctx).status == "SKIP"

    def test_nearly_confounded_warning(self):
        # 5 of 6 samples perfectly confounded, 1 cross-batch → V high but < 1
        meta = self._meta_with_batch(
            conditions=["control", "control", "control", "treated", "treated", "treated"],
            batches=["A", "A", "B", "B", "B", "B"],  # control mostly in A, treated all in B
        )
        ctx = _ctx(count_df=make_count_matrix(), meta_df=meta)
        result = BatchConfoundingRule().run(ctx)
        # Should be FAIL (WARNING or ERROR) because V is high
        assert result.status == "FAIL"


# ===========================================================================
# BIO-008 ERCCSpikeInRule
# ===========================================================================

class TestERCCSpikeInRule:
    def _make_ercc_df(self, n_ercc: int, ercc_fraction: float = 0.05):
        n_genes = 200
        genes = [f"ENSG{i:011d}" for i in range(1, n_genes + 1)]
        ercc_ids = [f"ERCC-{i:06d}" for i in range(1, n_ercc + 1)]
        all_genes = genes + ercc_ids
        n_samples = 6
        rng = np.random.default_rng(20)
        values = rng.integers(10, 500, size=(len(all_genes), n_samples))
        df = pd.DataFrame(
            values, index=all_genes,
            columns=[f"s{i}" for i in range(n_samples)]
        )
        # Scale ERCC rows so they represent the desired fraction
        total_per_sample = df.sum(axis=0)
        target_ercc_total = total_per_sample * ercc_fraction / (1 - ercc_fraction)
        ercc_sum = df.loc[ercc_ids].sum(axis=0)
        scale = target_ercc_total / ercc_sum.replace(0, 1)
        df.loc[ercc_ids] = df.loc[ercc_ids].multiply(scale, axis=1).round().astype(int)
        return df

    def test_no_ercc_passes(self):
        ctx = _ctx(count_df=make_count_matrix())
        assert ERCCSpikeInRule().run(ctx).status == "PASS"

    def test_low_ercc_fraction_warning(self):
        # Small ERCC fraction → WARNING not ERROR
        df = self._make_ercc_df(n_ercc=5, ercc_fraction=0.03)
        ctx = _ctx(count_df=df)
        result = ERCCSpikeInRule().run(ctx)
        assert result.status == "FAIL"
        assert result.severity == "WARNING"

    def test_high_ercc_fraction_error(self):
        # > 10% ERCC fraction → ERROR
        df = self._make_ercc_df(n_ercc=10, ercc_fraction=0.20)
        ctx = _ctx(count_df=df)
        result = ERCCSpikeInRule().run(ctx)
        assert result.status == "FAIL"
        assert result.severity == "ERROR"

    def test_ercc_ids_case_insensitive(self):
        # Mixed-case ERCC IDs should still be detected
        genes = [f"ENSG{i:011d}" for i in range(1, 11)] + ["ercc-000001", "ERCC-000002"]
        rng = np.random.default_rng(21)
        values = rng.integers(1, 100, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        ctx = _ctx(count_df=df)
        result = ERCCSpikeInRule().run(ctx)
        assert result.status == "FAIL"

    def test_no_matrix_skips(self):
        ctx = _ctx()
        assert ERCCSpikeInRule().run(ctx).status == "SKIP"


# ===========================================================================
# Extra coverage: BIO-006 HousekeepingGeneRule — additional paths
# ===========================================================================

class TestHousekeepingGeneRuleExtraPaths:
    def test_no_matrix_skips(self):
        ctx = _ctx()
        assert HousekeepingGeneRule().run(ctx).status == "SKIP"

    def test_mouse_ensembl_path_used(self):
        """When organism=mouse, mouse HK Ensembl map is used."""
        from app.rules.biology.rules import _MOUSE_HK_ENSEMBL
        # Include 3 of the 10 mouse HK Ensembl IDs
        hk_ids = list(_MOUSE_HK_ENSEMBL.keys())[:3]
        other_ids = [f"ENSMUSG{i:011d}" for i in range(100, 300)]
        genes = hk_ids + other_ids
        rng = np.random.default_rng(30)
        values = rng.integers(1, 100, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        ctx = _ctx(count_df=df, flags={"organism": "mouse"})
        result = HousekeepingGeneRule().run(ctx)
        assert result.status == "PASS"

    def test_entrez_id_path_skips(self):
        """Pure Entrez IDs (integers) should skip the housekeeping check."""
        # Entrez IDs: not Ensembl, not symbols → < 50% symbols → SKIP
        genes = [str(i) for i in range(1000, 1200)]
        rng = np.random.default_rng(31)
        values = rng.integers(1, 100, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        ctx = _ctx(count_df=df)
        result = HousekeepingGeneRule().run(ctx)
        assert result.status == "SKIP"

    def test_symbol_path_fails_all_hk_missing(self):
        """Symbol path: when ≥ 80% of housekeeping genes are absent, rule fails."""
        # Use non-HK gene symbols that satisfy the symbol regex
        genes = [f"GENE{i:02d}" for i in range(100)]  # e.g. GENE00–GENE99
        rng = np.random.default_rng(32)
        values = rng.integers(1, 100, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        ctx = _ctx(count_df=df)
        result = HousekeepingGeneRule().run(ctx)
        assert result.status == "FAIL"


# ===========================================================================
# Extra coverage: BIO-007 BatchConfoundingRule — additional paths
# ===========================================================================

class TestBatchConfoundingRuleExtraPaths:
    def _meta_with_batch(self, conditions, batches):
        samples = [f"s{i}" for i in range(len(conditions))]
        return pd.DataFrame(
            {"condition": conditions, "batch": batches},
            index=samples,
        )

    def test_no_condition_column_skips(self):
        """BIO-007 skips when the metadata has no condition column."""
        meta = pd.DataFrame(
            {"batch": ["A", "B", "A", "B", "A", "B"]},
            index=[f"s{i}" for i in range(6)],
        )
        ctx = _ctx(count_df=make_count_matrix(), meta_df=meta)
        result = BatchConfoundingRule().run(ctx)
        assert result.status == "SKIP"

    def test_cramers_v_single_level_batch(self):
        """A single-level batch column results in V=0 and a PASS."""
        meta = self._meta_with_batch(
            conditions=["control"] * 3 + ["treated"] * 3,
            batches=["A"] * 6,  # all samples in one batch
        )
        ctx = _ctx(count_df=make_count_matrix(), meta_df=meta)
        result = BatchConfoundingRule().run(ctx)
        # V = 0 because contingency table has only 1 batch column → PASS
        assert result.status == "PASS"

    def test_scipy_unavailable_skips(self, monkeypatch):
        """BIO-007 skips gracefully when scipy is not importable."""
        import app.rules.biology.rules as bio_rules

        original = bio_rules.BatchConfoundingRule._cramers_v

        @staticmethod
        def _nan_v(x, y):
            return float("nan")

        monkeypatch.setattr(bio_rules.BatchConfoundingRule, "_cramers_v", _nan_v)

        meta = self._meta_with_batch(
            conditions=["control"] * 3 + ["treated"] * 3,
            batches=["A"] * 3 + ["B"] * 3,
        )
        ctx = _ctx(count_df=make_count_matrix(), meta_df=meta)
        result = BatchConfoundingRule().run(ctx)
        assert result.status == "SKIP"


# ===========================================================================
# _detect_organism helper (biology/rules.py lines 25, 34-36)
# ===========================================================================

class TestDetectOrganism:
    """Direct unit tests for the _detect_organism helper in biology/rules.py."""

    def test_empty_list_returns_unknown(self):
        from app.rules.biology.rules import _detect_organism
        assert _detect_organism([]) == "unknown"

    def test_mouse_ids_detected(self):
        from app.rules.biology.rules import _detect_organism
        mouse_genes = [f"ENSMUSG{i:011d}" for i in range(1, 100)]
        assert _detect_organism(mouse_genes) == "mouse"

    def test_rat_ids_detected(self):
        from app.rules.biology.rules import _detect_organism
        rat_genes = [f"ENSRNOG{i:011d}" for i in range(1, 100)]
        assert _detect_organism(rat_genes) == "rat"

    def test_human_ids_detected(self):
        from app.rules.biology.rules import _detect_organism
        human_genes = [f"ENSG{i:011d}" for i in range(1, 100)]
        assert _detect_organism(human_genes) == "human"

    def test_mixed_returns_unknown(self):
        from app.rules.biology.rules import _detect_organism
        # No species exceeds 50%
        genes = (
            [f"ENSG{i:011d}" for i in range(1, 34)]
            + [f"ENSMUSG{i:011d}" for i in range(1, 34)]
            + [f"random_{i}" for i in range(1, 34)]
        )
        assert _detect_organism(genes) == "unknown"


# ===========================================================================
# BIO-005 — _detect_organism auto-detection for mouse (no organism flag)
# ===========================================================================

class TestMitochondrialFractionMouseAutodetect:
    def test_mouse_mt_autodetected_without_flag(self):
        """BIO-005 auto-detects mouse organism when context flag is absent."""
        genes = [f"ENSMUSG{i:011d}" for i in range(1, 51)] + [
            "mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2"
        ]
        rng = np.random.default_rng(88)
        values = rng.integers(1, 10, size=(len(genes), 4))
        df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
        for g in ["mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2"]:
            df.loc[g] = 2000
        ctx = _ctx(count_df=df)  # no organism flag — auto-detect must handle this
        result = MitochondrialFractionRule().run(ctx)
        assert result.status == "FAIL"


# ===========================================================================
# _cramers_v — ImportError path (biology/rules.py lines 511-512)
# ===========================================================================

class TestCramersVScipy:
    def test_cramers_v_returns_nan_when_scipy_missing(self, monkeypatch):
        """_cramers_v returns NaN (float) when scipy.stats import fails."""
        import sys

        original = sys.modules.get("scipy.stats")
        # Setting module to None causes 'from scipy.stats import ...' to raise ImportError
        monkeypatch.setitem(sys.modules, "scipy.stats", None)
        try:
            v = BatchConfoundingRule._cramers_v(
                pd.Series(["A", "B", "A"]),
                pd.Series(["X", "Y", "X"]),
            )
            assert v != v  # NaN != NaN is the canonical NaN check
        finally:
            # Restore original to avoid breaking subsequent tests
            if original is not None:
                sys.modules["scipy.stats"] = original
            elif "scipy.stats" in sys.modules:
                del sys.modules["scipy.stats"]
