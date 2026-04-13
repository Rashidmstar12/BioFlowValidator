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
