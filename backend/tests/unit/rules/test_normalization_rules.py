"""Unit tests for normalization validation rules."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from app.rules.normalization.rules import (
    AllZeroGeneRule,
    DuplicateRowRule,
    LibrarySizeRule,
    LowCountDominanceRule,
    NonIntegerCountRule,
    ZeroLibraryRule,
)
from tests.conftest import _ctx_from_dfs, make_count_matrix


# ── NonIntegerCountRule ───────────────────────────────────────────────────────

def test_non_integer_fails_float_counts():
    import numpy as np
    rng = np.random.default_rng(7)
    # All values are non-integer floats
    values = rng.uniform(0.1, 999.9, size=(20, 6))
    df = pd.DataFrame(values, index=[f"G{i}" for i in range(20)], columns=[f"s{i}" for i in range(6)])
    ctx = _ctx_from_dfs(count_df=df)
    assert NonIntegerCountRule().run(ctx).status == "FAIL"


def test_non_integer_passes_int_counts():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert NonIntegerCountRule().run(ctx).status == "PASS"


def test_non_integer_fails_tpm_like():
    # TPM-like values (non-integer, 0-10000 range)
    rng = np.random.default_rng(0)
    values = rng.uniform(0, 1000, size=(20, 6))
    df = pd.DataFrame(values, index=[f"G{i}" for i in range(20)], columns=[f"s{i}" for i in range(6)])
    ctx = _ctx_from_dfs(count_df=df)
    assert NonIntegerCountRule().run(ctx).status == "FAIL"


# ── LibrarySizeRule ───────────────────────────────────────────────────────────

def test_library_size_fails_extreme_ratio():
    rng = np.random.default_rng(1)
    df = make_count_matrix()
    df["ctrl_1"] = df["ctrl_1"] * 100  # 100x more counts
    ctx = _ctx_from_dfs(count_df=df)
    result = LibrarySizeRule().run(ctx)
    assert result.status == "FAIL"


def test_library_size_passes_balanced():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    result = LibrarySizeRule().run(ctx)
    assert result.status == "PASS"


# ── ZeroLibraryRule ───────────────────────────────────────────────────────────

def test_zero_library_fails():
    df = make_count_matrix()
    df["ctrl_1"] = 0
    ctx = _ctx_from_dfs(count_df=df)
    assert ZeroLibraryRule().run(ctx).status == "FAIL"


def test_zero_library_passes():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert ZeroLibraryRule().run(ctx).status == "PASS"


# ── AllZeroGeneRule ───────────────────────────────────────────────────────────

def test_all_zero_gene_fails():
    df = make_count_matrix()
    df.iloc[0] = 0
    df.iloc[1] = 0
    ctx = _ctx_from_dfs(count_df=df)
    assert AllZeroGeneRule().run(ctx).status == "FAIL"


def test_all_zero_gene_passes():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert AllZeroGeneRule().run(ctx).status == "PASS"


# ── LowCountDominanceRule ─────────────────────────────────────────────────────

def test_low_count_dominance_fails():
    # 90% of genes have median < 1
    df = pd.DataFrame(
        np.zeros((20, 6)),
        index=[f"G{i}" for i in range(20)],
        columns=[f"s{i}" for i in range(6)],
    )
    df.iloc[:2] = 100  # only 2 genes have counts
    ctx = _ctx_from_dfs(count_df=df)
    assert LowCountDominanceRule().run(ctx).status == "FAIL"


def test_low_count_dominance_passes():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert LowCountDominanceRule().run(ctx).status == "PASS"


# ── DuplicateRowRule ──────────────────────────────────────────────────────────

def test_duplicate_row_fails():
    df = make_count_matrix()
    df.iloc[1] = df.iloc[0]  # row 1 = row 0
    ctx = _ctx_from_dfs(count_df=df)
    assert DuplicateRowRule().run(ctx).status == "FAIL"


def test_duplicate_row_passes():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert DuplicateRowRule().run(ctx).status == "PASS"


# ── _skip when count_matrix is None ──────────────────────────────────────────

def test_non_integer_skips_no_matrix():
    ctx = _ctx_from_dfs()
    assert NonIntegerCountRule().run(ctx).status == "SKIP"


def test_non_integer_skips_all_nan():
    """A matrix that coerces entirely to NaN has no numeric values to check."""
    df = pd.DataFrame(
        {"s1": ["abc", "def"], "s2": ["xyz", "uvw"]},
        index=["G1", "G2"],
    )
    ctx = _ctx_from_dfs(count_df=df)
    assert NonIntegerCountRule().run(ctx).status == "SKIP"


def test_library_size_skips_no_matrix():
    ctx = _ctx_from_dfs()
    assert LibrarySizeRule().run(ctx).status == "SKIP"


def test_library_size_skips_too_few_positive():
    """When all samples have near-zero libraries and only 1 positive remains, skip."""
    df = pd.DataFrame(
        {"s1": [0, 0, 0], "s2": [0, 0, 0]},
        index=["G1", "G2", "G3"],
    )
    ctx = _ctx_from_dfs(count_df=df)
    # lib_sizes.min() == 0 → goes into the positive-filter branch
    # both samples have 0 counts → len(positive) < 2 → SKIP
    result = LibrarySizeRule().run(ctx)
    assert result.status == "SKIP"


def test_library_size_with_one_zero_sample_fails_ratio():
    """LibrarySizeRule computes ratio from positive-only samples when one is zero."""
    df = make_count_matrix()
    # Set one sample to all zeros, keep a 100× imbalance among the rest
    df["ctrl_1"] = 0
    df["ctrl_2"] = df["ctrl_2"] * 200  # creates >10× ratio among positives
    ctx = _ctx_from_dfs(count_df=df)
    result = LibrarySizeRule().run(ctx)
    # With one zero sample, uses ratio from the 5 positives; 200× > 10×
    assert result.status == "FAIL"


def test_zero_library_skips_no_matrix():
    ctx = _ctx_from_dfs()
    assert ZeroLibraryRule().run(ctx).status == "SKIP"


def test_all_zero_gene_skips_no_matrix():
    ctx = _ctx_from_dfs()
    assert AllZeroGeneRule().run(ctx).status == "SKIP"


def test_low_count_dominance_skips_no_matrix():
    ctx = _ctx_from_dfs()
    assert LowCountDominanceRule().run(ctx).status == "SKIP"


def test_duplicate_row_skips_no_matrix():
    ctx = _ctx_from_dfs()
    assert DuplicateRowRule().run(ctx).status == "SKIP"
