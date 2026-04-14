"""Unit tests for biology sanity check rules."""
from __future__ import annotations

import pandas as pd
import pytest

from app.rules.biology.rules import (
    ConditionLabelSanityRule,
    HighCountGeneRule,
    MetadataCardinalityRule,
    MitochondrialFractionRule,
    SingleConditionRule,
)
from tests.conftest import _ctx_from_dfs, make_count_matrix, make_metadata


# ── SingleConditionRule ───────────────────────────────────────────────────────

def test_single_condition_fails():
    meta = make_metadata(
        samples=["s1", "s2", "s3"],
        conditions=["control", "control", "control"],
    )
    ctx = _ctx_from_dfs(count_df=make_count_matrix(samples=["s1", "s2", "s3"]), meta_df=meta)
    result = SingleConditionRule().run(ctx)
    assert result.status == "FAIL"
    assert result.severity == "ERROR"


def test_single_condition_passes():
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=make_metadata())
    assert SingleConditionRule().run(ctx).status == "PASS"


def test_single_condition_skips_no_metadata():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert SingleConditionRule().run(ctx).status == "SKIP"


# ── ConditionLabelSanityRule ──────────────────────────────────────────────────

def test_numeric_labels_fail():
    meta = pd.DataFrame(
        {"condition": ["1", "1", "1", "2", "2", "2"]},
        index=["s1", "s2", "s3", "s4", "s5", "s6"],
    )
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=meta)
    result = ConditionLabelSanityRule().run(ctx)
    assert result.status == "FAIL"


def test_string_labels_pass():
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=make_metadata())
    assert ConditionLabelSanityRule().run(ctx).status == "PASS"


# ── MetadataCardinalityRule ───────────────────────────────────────────────────

def test_high_cardinality_fails():
    meta = pd.DataFrame(
        {"condition": [f"sample_{i}" for i in range(6)]},
        index=[f"s{i}" for i in range(6)],
    )
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=meta)
    result = MetadataCardinalityRule().run(ctx)
    assert result.status == "FAIL"


def test_normal_cardinality_passes():
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=make_metadata())
    assert MetadataCardinalityRule().run(ctx).status == "PASS"


# ── HighCountGeneRule ─────────────────────────────────────────────────────────

def test_high_count_gene_fails():
    import numpy as np
    # One gene dominates > 50% of each sample's library
    df = pd.DataFrame(
        np.ones((20, 6), dtype=int) * 10,
        index=[f"G{i}" for i in range(20)],
        columns=[f"s{i}" for i in range(6)],
    )
    df.iloc[0] = 10_000  # this gene = ~98% of library (other 19 genes = 10 each = 190 total)
    ctx = _ctx_from_dfs(count_df=df)
    result = HighCountGeneRule().run(ctx)
    assert result.status == "FAIL"


def test_high_count_gene_passes():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert HighCountGeneRule().run(ctx).status == "PASS"


# ── MitochondrialFractionRule ─────────────────────────────────────────────────

def test_mt_fraction_fails():
    import numpy as np
    genes = [f"ENSG{i:011d}" for i in range(1, 11)] + ["MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2"]
    rng = np.random.default_rng(5)
    values = rng.integers(1, 10, size=(len(genes), 6))
    df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(6)])
    # Make MT genes dominate
    df.loc["MT-ND1"] = 1000
    df.loc["MT-ND2"] = 1000
    df.loc["MT-CO1"] = 1000
    df.loc["MT-CO2"] = 1000
    ctx = _ctx_from_dfs(count_df=df)
    result = MitochondrialFractionRule().run(ctx)
    assert result.status == "FAIL"


def test_mt_fraction_passes_no_mt_genes():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    # ENSG IDs — no MT pattern, should PASS
    result = MitochondrialFractionRule().run(ctx)
    assert result.status == "PASS"


# ── _skip / edge cases for biology rules ─────────────────────────────────────

def test_single_condition_uses_fallback_column():
    """BIO-001 must use first non-numeric column when no standard name matches."""
    idx = [f"s{i}" for i in range(6)]
    meta = pd.DataFrame(
        {"tissue": pd.Series(["liver", "liver", "liver", "kidney", "kidney", "kidney"], dtype=object, index=idx)},
        index=idx,
    )
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=meta)
    result = SingleConditionRule().run(ctx)
    assert result.status == "PASS"


def test_single_condition_skips_no_condition_column():
    """BIO-001 must skip when no non-numeric column exists."""
    meta = pd.DataFrame(
        {"age": [30, 40, 50, 60, 70, 80]},
        index=[f"s{i}" for i in range(6)],
    )
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=meta)
    assert SingleConditionRule().run(ctx).status == "SKIP"


def test_condition_label_sanity_skips_no_metadata():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert ConditionLabelSanityRule().run(ctx).status == "SKIP"


def test_condition_label_sanity_skips_no_condition_column():
    meta = pd.DataFrame({"age": [30, 40, 50]}, index=["s1", "s2", "s3"])
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=meta)
    assert ConditionLabelSanityRule().run(ctx).status == "SKIP"


def test_condition_label_sanity_uses_fallback_column():
    """BIO-002 must fall back to the first non-numeric column."""
    idx = [f"s{i}" for i in range(6)]
    meta = pd.DataFrame(
        {"tissue": pd.Series(["liver", "liver", "liver", "kidney", "kidney", "kidney"], dtype=object, index=idx)},
        index=idx,
    )
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=meta)
    result = ConditionLabelSanityRule().run(ctx)
    assert result.status == "PASS"


def test_metadata_cardinality_skips_no_metadata():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert MetadataCardinalityRule().run(ctx).status == "SKIP"


def test_metadata_cardinality_skips_few_samples():
    """BIO-003 must skip when fewer than 3 samples are present."""
    meta = pd.DataFrame(
        {"condition": ["control", "treated"]},
        index=["s1", "s2"],
    )
    ctx = _ctx_from_dfs(count_df=make_count_matrix(samples=["s1", "s2"]), meta_df=meta)
    assert MetadataCardinalityRule().run(ctx).status == "SKIP"


def test_high_count_gene_skips_no_matrix():
    ctx = _ctx_from_dfs()
    assert HighCountGeneRule().run(ctx).status == "SKIP"


def test_high_count_gene_skips_all_zero_library():
    """BIO-004 must skip when every sample has zero total counts."""
    import numpy as np
    df = pd.DataFrame(
        np.zeros((10, 4), dtype=int),
        index=[f"G{i}" for i in range(10)],
        columns=[f"s{i}" for i in range(4)],
    )
    ctx = _ctx_from_dfs(count_df=df)
    assert HighCountGeneRule().run(ctx).status == "SKIP"


def test_mt_fraction_skips_no_matrix():
    ctx = _ctx_from_dfs()
    assert MitochondrialFractionRule().run(ctx).status == "SKIP"


def test_mt_fraction_passes_below_threshold():
    """BIO-005 passes when MT genes are present but fraction is below 30%."""
    import numpy as np
    genes = [f"ENSG{i:011d}" for i in range(1, 97)] + ["MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2"]
    rng = np.random.default_rng(99)
    values = rng.integers(100, 1000, size=(len(genes), 4))
    df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(4)])
    # Keep MT genes at very low counts (< 5% of total)
    for g in ["MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2"]:
        df.loc[g] = 5
    ctx = _ctx_from_dfs(count_df=df)
    ctx.flags["organism"] = "human"
    result = MitochondrialFractionRule().run(ctx)
    assert result.status == "PASS"
