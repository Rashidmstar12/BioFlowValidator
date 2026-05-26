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


def test_non_integer_tpm_sums_to_one_million_is_error():
    genes = [f"G{i}" for i in range(1100)]
    values = np.random.uniform(0.1, 1000.0, size=(len(genes), 6))
    values = values / values.sum(axis=0) * 1_000_000
    df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(6)])
    ctx = _ctx_from_dfs(count_df=df)
    result = NonIntegerCountRule().run(ctx)
    assert result.status == "FAIL"
    assert result.severity == "ERROR"


def test_non_integer_expected_counts_is_warning():
    genes = [f"G{i}" for i in range(1100)]
    values = np.random.uniform(0.1, 10000.0, size=(len(genes), 6))
    df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(6)])
    ctx = _ctx_from_dfs(count_df=df)
    result = NonIntegerCountRule().run(ctx)
    assert result.status == "FAIL"
    assert result.severity == "WARNING"


def test_non_integer_fpkm_is_error():
    genes = [f"G{i}" for i in range(1100)]
    values = np.random.uniform(0.1, 10.0, size=(len(genes), 6))
    df = pd.DataFrame(values, index=genes, columns=[f"s{i}" for i in range(6)])
    ctx = _ctx_from_dfs(count_df=df)
    result = NonIntegerCountRule().run(ctx)
    assert result.status == "FAIL"
    assert result.severity == "ERROR"
