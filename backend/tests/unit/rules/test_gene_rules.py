"""Unit tests for gene ID validation rules."""
from __future__ import annotations

import pandas as pd
import pytest

from app.rules.gene.rules import (
    DuplicateGeneRule,
    GeneIDFormatRule,
    NonBiologicalIDRule,
    VersionSuffixRule,
)
from tests.conftest import _ctx_from_dfs, make_count_matrix


def _make_ctx(genes: list[str]) -> object:
    df = make_count_matrix(genes=genes)
    return _ctx_from_dfs(count_df=df)


# ── GeneIDFormatRule ──────────────────────────────────────────────────────────

def test_gene_id_format_passes_ensembl():
    genes = [f"ENSG{i:011d}" for i in range(1, 11)]
    assert GeneIDFormatRule().run(_make_ctx(genes)).status == "PASS"


def test_gene_id_format_fails_mixed():
    ensembl = [f"ENSG{i:011d}" for i in range(1, 6)]
    symbols = ["ACTB", "GAPDH", "TP53", "BRCA1", "MYC"]
    ctx = _make_ctx(ensembl + symbols)
    result = GeneIDFormatRule().run(ctx)
    assert result.status == "FAIL"
    assert "ensembl" in str(result.details)


def test_gene_id_format_passes_symbols():
    genes = ["ACTB", "GAPDH", "TP53", "BRCA1", "MYC", "EGFR", "CDK2", "RAF1", "MAP2K1", "PIK3CA"]
    assert GeneIDFormatRule().run(_make_ctx(genes)).status == "PASS"


# ── DuplicateGeneRule ─────────────────────────────────────────────────────────

def test_duplicate_gene_fails():
    genes = ["ENSG00000000001", "ENSG00000000001", "ENSG00000000002"]
    df = pd.DataFrame([[1, 2], [3, 4], [5, 6]], index=genes, columns=["s1", "s2"])
    ctx = _ctx_from_dfs(count_df=df)
    result = DuplicateGeneRule().run(ctx)
    assert result.status == "FAIL"


def test_duplicate_gene_passes():
    ctx = _make_ctx([f"ENSG{i:011d}" for i in range(1, 11)])
    assert DuplicateGeneRule().run(ctx).status == "PASS"


# ── VersionSuffixRule ─────────────────────────────────────────────────────────

def test_version_suffix_fails():
    genes = [f"ENSG{i:011d}.3" for i in range(1, 11)]
    ctx = _make_ctx(genes)
    result = VersionSuffixRule().run(ctx)
    assert result.status == "FAIL"


def test_version_suffix_passes():
    genes = [f"ENSG{i:011d}" for i in range(1, 11)]
    assert VersionSuffixRule().run(_make_ctx(genes)).status == "PASS"


# ── NonBiologicalIDRule ───────────────────────────────────────────────────────

def test_non_biological_fails():
    genes = [f"random_xyz_{i}" for i in range(10)]
    ctx = _make_ctx(genes)
    result = NonBiologicalIDRule().run(ctx)
    assert result.status == "FAIL"
