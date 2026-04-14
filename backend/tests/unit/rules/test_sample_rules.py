"""Unit tests for sample validation rules."""
from __future__ import annotations

import pandas as pd
import pytest

from app.rules.sample.rules import (
    DuplicateSampleRule,
    NearIdenticalSampleRule,
    ReplicateCountRule,
    SampleMatchRule,
    SampleOrderRule,
)
from tests.conftest import _ctx_from_dfs, make_count_matrix, make_metadata


# ── SampleMatchRule ───────────────────────────────────────────────────────────

def test_sample_match_passes():
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=make_metadata())
    assert SampleMatchRule().run(ctx).status == "PASS"


def test_sample_match_fails_extra_in_matrix():
    count = make_count_matrix(samples=["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3", "extra"])
    meta = make_metadata()
    ctx = _ctx_from_dfs(count_df=count, meta_df=meta)
    result = SampleMatchRule().run(ctx)
    assert result.status == "FAIL"
    assert "extra" in str(result.details)


def test_sample_match_fails_extra_in_meta():
    count = make_count_matrix()
    meta = make_metadata(
        samples=["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3", "ghost"],
        conditions=["control"] * 3 + ["treated"] * 3 + ["treated"],
    )
    ctx = _ctx_from_dfs(count_df=count, meta_df=meta)
    result = SampleMatchRule().run(ctx)
    assert result.status == "FAIL"


def test_sample_match_skips_no_metadata():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert SampleMatchRule().run(ctx).status == "SKIP"


# ── DuplicateSampleRule ───────────────────────────────────────────────────────

def test_duplicate_sample_fails_in_matrix():
    df = make_count_matrix()
    df.columns = ["ctrl_1", "ctrl_1", "ctrl_3", "treat_1", "treat_2", "treat_3"]
    ctx = _ctx_from_dfs(count_df=df)
    assert DuplicateSampleRule().run(ctx).status == "FAIL"


def test_duplicate_sample_passes():
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=make_metadata())
    assert DuplicateSampleRule().run(ctx).status == "PASS"


# ── ReplicateCountRule ────────────────────────────────────────────────────────

def test_replicate_count_error_one_replicate():
    meta = make_metadata(
        samples=["ctrl_1", "treat_1"],
        conditions=["control", "treated"],
    )
    ctx = _ctx_from_dfs(count_df=make_count_matrix(samples=["ctrl_1", "treat_1"], genes=["G1"]), meta_df=meta)
    result = ReplicateCountRule().run(ctx)
    assert result.status == "FAIL"
    assert result.severity == "ERROR"


def test_replicate_count_warning_two_replicates():
    meta = make_metadata(
        samples=["ctrl_1", "ctrl_2", "treat_1", "treat_2"],
        conditions=["control", "control", "treated", "treated"],
    )
    ctx = _ctx_from_dfs(
        count_df=make_count_matrix(samples=["ctrl_1", "ctrl_2", "treat_1", "treat_2"]),
        meta_df=meta,
    )
    result = ReplicateCountRule().run(ctx)
    assert result.status == "FAIL"
    assert result.severity == "WARNING"


def test_replicate_count_passes_three():
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=make_metadata())
    assert ReplicateCountRule().run(ctx).status == "PASS"


# ── NearIdenticalSampleRule ───────────────────────────────────────────────────

def test_near_identical_fails():
    import numpy as np
    df = make_count_matrix()
    df["ctrl_3"] = df["ctrl_1"]  # exact duplicate
    ctx = _ctx_from_dfs(count_df=df)
    result = NearIdenticalSampleRule().run(ctx)
    assert result.status == "FAIL"


def test_near_identical_passes():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert NearIdenticalSampleRule().run(ctx).status == "PASS"


# ── _skip edge cases ──────────────────────────────────────────────────────────

def test_sample_match_skips_no_matrix():
    """SMP-001 must skip when count matrix is absent."""
    ctx = _ctx_from_dfs(meta_df=make_metadata())
    assert SampleMatchRule().run(ctx).status == "SKIP"


def test_sample_order_skips_no_matrix_or_meta():
    """SMP-002 must skip when either file is absent."""
    ctx = _ctx_from_dfs()
    assert SampleOrderRule().run(ctx).status == "SKIP"


def test_sample_order_passes_when_order_differs_but_only_subset_common():
    """SMP-002 checks only common samples; all-common but different order → FAIL."""
    count = make_count_matrix()
    # Reorder metadata rows compared to count matrix columns
    meta = make_metadata(
        samples=["treat_1", "treat_2", "treat_3", "ctrl_1", "ctrl_2", "ctrl_3"],
        conditions=["treated"] * 3 + ["control"] * 3,
    )
    ctx = _ctx_from_dfs(count_df=count, meta_df=meta)
    result = SampleOrderRule().run(ctx)
    assert result.status == "FAIL"


def test_duplicate_sample_fails_in_metadata():
    """SMP-003 must catch duplicate row indices in metadata."""
    meta = pd.DataFrame(
        {"condition": ["control", "control", "treated"]},
        index=["s1", "s1", "s2"],
    )
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=meta)
    result = DuplicateSampleRule().run(ctx)
    assert result.status == "FAIL"


def test_replicate_count_skips_no_metadata():
    ctx = _ctx_from_dfs(count_df=make_count_matrix())
    assert ReplicateCountRule().run(ctx).status == "SKIP"


def test_replicate_count_skips_no_condition_column():
    """SMP-004 must skip when no condition-like column exists."""
    meta = pd.DataFrame(
        {"age": [30, 40, 50, 60, 70, 80]},
        index=[f"s{i}" for i in range(6)],
    )
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=meta)
    assert ReplicateCountRule().run(ctx).status == "SKIP"


def test_replicate_count_uses_fallback_condition_column():
    """SMP-004 falls back to first non-numeric column when no standard name matches."""
    idx = [f"s{i}" for i in range(6)]
    meta = pd.DataFrame(
        {"tissue": pd.Series(["liver"] * 3 + ["kidney"] * 3, dtype=object, index=idx)},
        index=idx,
    )
    ctx = _ctx_from_dfs(count_df=make_count_matrix(), meta_df=meta)
    result = ReplicateCountRule().run(ctx)
    assert result.status == "PASS"


def test_near_identical_skips_no_matrix():
    ctx = _ctx_from_dfs()
    assert NearIdenticalSampleRule().run(ctx).status == "SKIP"


def test_near_identical_skips_single_sample():
    """SMP-005 must skip when the count matrix has only 1 sample."""
    df = make_count_matrix(samples=["only_sample"])
    ctx = _ctx_from_dfs(count_df=df)
    assert NearIdenticalSampleRule().run(ctx).status == "SKIP"


def test_near_identical_skips_over_500_samples():
    """SMP-005 must skip when sample count exceeds 500 to avoid O(n²) cost."""
    import numpy as np
    n_samples = 501
    df = pd.DataFrame(
        np.ones((10, n_samples), dtype=int),
        index=[f"G{i}" for i in range(10)],
        columns=[f"s{i}" for i in range(n_samples)],
    )
    ctx = _ctx_from_dfs(count_df=df)
    result = NearIdenticalSampleRule().run(ctx)
    assert result.status == "SKIP"
