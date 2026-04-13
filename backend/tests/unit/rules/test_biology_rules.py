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
