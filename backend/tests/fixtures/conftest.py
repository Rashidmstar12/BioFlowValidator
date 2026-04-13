"""Shared pytest fixtures for unit and integration tests."""
from __future__ import annotations

import io
import textwrap

import pandas as pd
import pytest

from app.models.context import ValidationContext


# ── helpers ──────────────────────────────────────────────────────────────────

def _ctx_from_dfs(
    count_df: pd.DataFrame | None = None,
    meta_df: pd.DataFrame | None = None,
    count_bytes: bytes = b"",
    meta_bytes: bytes = b"",
) -> ValidationContext:
    ctx = ValidationContext(
        count_matrix_bytes=count_bytes,
        metadata_bytes=meta_bytes,
        count_matrix=count_df,
        metadata=meta_df,
        count_delimiter="\t",
        metadata_delimiter="\t",
        count_filename="counts.tsv",
        metadata_filename="metadata.tsv",
    )
    return ctx


def make_count_matrix(
    genes: list[str] | None = None,
    samples: list[str] | None = None,
    values=None,
) -> pd.DataFrame:
    if genes is None:
        genes = [f"ENSG{i:011d}" for i in range(1, 21)]
    if samples is None:
        samples = ["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"]
    if values is None:
        import numpy as np
        rng = np.random.default_rng(42)
        values = rng.integers(0, 1000, size=(len(genes), len(samples)))
    df = pd.DataFrame(values, index=genes, columns=samples)
    df.index.name = "gene_id"
    return df


def make_metadata(
    samples: list[str] | None = None,
    conditions: list[str] | None = None,
) -> pd.DataFrame:
    if samples is None:
        samples = ["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"]
    if conditions is None:
        conditions = ["control"] * 3 + ["treated"] * 3
    df = pd.DataFrame({"condition": conditions}, index=samples)
    df.index.name = "sample_id"
    return df


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def valid_count_matrix() -> pd.DataFrame:
    return make_count_matrix()


@pytest.fixture
def valid_metadata() -> pd.DataFrame:
    return make_metadata()


@pytest.fixture
def valid_context(valid_count_matrix, valid_metadata) -> ValidationContext:
    count_tsv = valid_count_matrix.to_csv(sep="\t").encode()
    meta_tsv = valid_metadata.to_csv(sep="\t").encode()
    return _ctx_from_dfs(
        count_df=valid_count_matrix,
        meta_df=valid_metadata,
        count_bytes=count_tsv,
        meta_bytes=meta_tsv,
    )
