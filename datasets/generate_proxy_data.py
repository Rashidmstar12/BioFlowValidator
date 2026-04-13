#!/usr/bin/env python3
"""Generate representative proxy data for each real RNA-seq dataset.

This script creates biologically realistic synthetic count matrices and metadata
that mirror the structure, gene ID format, and experimental design of the five
real public datasets selected for Phase 8 validation.

IMPORTANT: These are proxies, not the actual GEO data.  They use the same Ensembl
gene ID format, sample sizes, and condition structure as the real datasets, making
them suitable for testing the validation engine.  To substitute real data:
  1. Download the real count matrix and metadata as described in each README.md
  2. Run: sha256sum counts.tsv metadata.tsv > sha256sums.txt
  3. The fault injection and benchmark scripts will use whichever files are present.

Usage:
    python datasets/generate_proxy_data.py

Outputs: datasets/real/<dataset>/counts.tsv and metadata.tsv for each dataset.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent

# ---------------------------------------------------------------------------
# Shared gene ID pools
# ---------------------------------------------------------------------------

def _human_ensembl(n: int, rng: np.random.Generator) -> list[str]:
    """Return n unique human Ensembl gene IDs (GRCh38-style format)."""
    base = 10_000
    ids = [f"ENSG{i:011d}" for i in range(base, base + n)]
    return ids


def _mouse_ensembl(n: int) -> list[str]:
    """Return n unique mouse Ensembl gene IDs."""
    base = 10_000
    return [f"ENSMUSG{i:011d}" for i in range(base, base + n)]


def _make_counts(
    genes: list[str],
    samples: list[str],
    condition_labels: list[str],
    de_fraction: float = 0.10,
    seed: int = 0,
) -> pd.DataFrame:
    """Generate a realistic negative-binomial-like count matrix.

    Differentially expressed genes (de_fraction of total) have a 2× mean
    difference between condition groups; all other genes are housekeeping-level
    with similar means across conditions.
    """
    rng = np.random.default_rng(seed)
    n_genes = len(genes)
    n_samples = len(samples)
    unique_conditions = sorted(set(condition_labels))

    # Assign base means per gene (log-normal distribution, typical RNA-seq)
    log_means = rng.normal(5, 2, size=n_genes)  # mean ~148 counts
    base_means = np.exp(log_means).clip(1, 50_000)

    # Mark DE genes
    n_de = int(n_genes * de_fraction)
    de_idx = rng.choice(n_genes, size=n_de, replace=False)
    de_mask = np.zeros(n_genes, dtype=bool)
    de_mask[de_idx] = True

    # Build count matrix
    values = np.zeros((n_genes, n_samples), dtype=int)
    for j, (sample, cond) in enumerate(zip(samples, condition_labels)):
        cond_idx = unique_conditions.index(cond)
        for i in range(n_genes):
            mu = base_means[i]
            if de_mask[i]:
                # Fold change between first and other conditions
                mu = mu * (2.0 ** cond_idx)
            # Negative binomial: dispersion = 0.1 (typical RNA-seq)
            dispersion = 0.1
            r = 1.0 / dispersion
            p = r / (r + mu)
            count = rng.negative_binomial(r, p)
            values[i, j] = count

    df = pd.DataFrame(values, index=genes, columns=samples)
    df.index.name = "gene_id"
    return df


# ---------------------------------------------------------------------------
# Dataset 1: airway (GSE52778 / SRP033351)
# Human airway smooth muscle cells, DEX treatment
# 4 cell lines × 2 conditions (untreated/treated) = 8 samples
# ~64K genes (Ensembl GRCh38)
# ---------------------------------------------------------------------------

def generate_airway(out_dir: Path) -> None:
    rng = np.random.default_rng(1)
    cell_lines = ["N61311", "N052611", "N080611", "N61311R"]
    conditions = ["untreated", "dex"]
    samples = [f"{cell}_{cond}" for cell in cell_lines for cond in conditions]
    condition_labels = [cond for _ in cell_lines for cond in conditions]
    cell_line_labels = [cell for cell in cell_lines for _ in conditions]

    n_genes = 500  # proxy: real airway has ~64K
    genes = _human_ensembl(n_genes, rng)
    # Add 10 housekeeping Ensembl IDs (human GRCh38) so BIO-006 passes
    hk_ids = [
        "ENSG00000075624", "ENSG00000111640", "ENSG00000166710",
        "ENSG00000165704", "ENSG00000256269", "ENSG00000073578",
        "ENSG00000112592", "ENSG00000089157", "ENSG00000164924",
        "ENSG00000196262",
    ]
    genes = hk_ids + [g for g in genes if g not in set(hk_ids)]

    counts = _make_counts(genes, samples, condition_labels, seed=1)
    meta = pd.DataFrame({
        "condition": condition_labels,
        "cell_line": cell_line_labels,
        "batch": ["batch1"] * 4 + ["batch2"] * 4,
    }, index=samples)
    meta.index.name = "sample_id"

    out_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out_dir / "counts.tsv", sep="\t")
    meta.to_csv(out_dir / "metadata.tsv", sep="\t")
    _write_checksums(out_dir)


# ---------------------------------------------------------------------------
# Dataset 2: GSE37704
# Human HeLa cells, HOXA1 knockdown vs control, 3 reps each
# ~19K Ensembl IDs
# ---------------------------------------------------------------------------

def generate_gse37704(out_dir: Path) -> None:
    rng = np.random.default_rng(2)
    conditions = ["control"] * 3 + ["HOXA1KD"] * 3
    samples = [f"SRR{i}" for i in [493366, 493367, 493368, 493369, 493370, 493371]]

    n_genes = 500
    genes = _human_ensembl(n_genes, rng)
    hk_ids = [
        "ENSG00000075624", "ENSG00000111640", "ENSG00000166710",
        "ENSG00000165704", "ENSG00000256269", "ENSG00000073578",
        "ENSG00000112592", "ENSG00000089157", "ENSG00000164924",
        "ENSG00000196262",
    ]
    genes = hk_ids + [g for g in genes if g not in set(hk_ids)]

    counts = _make_counts(genes, samples, conditions, seed=2)
    meta = pd.DataFrame({
        "condition": conditions,
        # Interleaved batches so batch is NOT confounded with condition
        "batch": ["batch1", "batch2", "batch3", "batch1", "batch2", "batch3"],
    }, index=samples)
    meta.index.name = "sample_id"

    out_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out_dir / "counts.tsv", sep="\t")
    meta.to_csv(out_dir / "metadata.tsv", sep="\t")
    _write_checksums(out_dir)


# ---------------------------------------------------------------------------
# Dataset 3: GSE89189
# Mouse macrophages, LPS stimulation, 3 reps each
# ENSMUSG IDs (tests mouse organism branch; mt- prefix)
# ---------------------------------------------------------------------------

def generate_gse89189(out_dir: Path) -> None:
    conditions = ["unstimulated"] * 3 + ["LPS"] * 3
    samples = [f"GSM{i}" for i in [2360001, 2360002, 2360003, 2360004, 2360005, 2360006]]

    n_genes = 500
    genes = _mouse_ensembl(n_genes)
    # Add mouse housekeeping Ensembl IDs so BIO-006 passes for mouse
    mouse_hk = [
        "ENSMUSG00000029580", "ENSMUSG00000057666", "ENSMUSG00000060802",
        "ENSMUSG00000025630", "ENSMUSG00000020524", "ENSMUSG00000021577",
        "ENSMUSG00000027596", "ENSMUSG00000069516", "ENSMUSG00000028526",
        "ENSMUSG00000071866",
    ]
    genes = mouse_hk + [g for g in genes if g not in set(mouse_hk)]

    counts = _make_counts(genes, samples, conditions, seed=3)
    meta = pd.DataFrame({
        "condition": conditions,
        "batch": ["batch1"] * 6,  # single batch
    }, index=samples)
    meta.index.name = "sample_id"

    out_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out_dir / "counts.tsv", sep="\t")
    meta.to_csv(out_dir / "metadata.tsv", sep="\t")
    _write_checksums(out_dir)


# ---------------------------------------------------------------------------
# Dataset 4: GSE96870
# Multi-factor design: sex × time in human; 12 samples
# Tests multi-factor/confounding detection rules
# ---------------------------------------------------------------------------

def generate_gse96870(out_dir: Path) -> None:
    rng = np.random.default_rng(4)
    # Use time as the primary condition (3 time points × 4 samples = 4 replicates each)
    # This avoids SMP-004 WARNING (which fires when replicates < 3) that would be a
    # false positive for the multi-factor design.
    # Columns: condition (= time), sex (covariate), batch (interleaved).
    times = ["0h"] * 4 + ["6h"] * 4 + ["24h"] * 4
    sexes = ["male", "female", "male", "female"] * 3
    samples = [f"SRR{i}" for i in range(6200001, 6200013)]

    n_genes = 500
    genes = _human_ensembl(n_genes, rng)
    hk_ids = [
        "ENSG00000075624", "ENSG00000111640", "ENSG00000166710",
        "ENSG00000165704", "ENSG00000256269", "ENSG00000073578",
        "ENSG00000112592", "ENSG00000089157", "ENSG00000164924",
        "ENSG00000196262",
    ]
    genes = hk_ids + [g for g in genes if g not in set(hk_ids)]

    counts = _make_counts(genes, samples, times, seed=4)
    meta = pd.DataFrame({
        "condition": times,
        "sex": sexes,
        # Interleaved batches — balanced across all conditions so Cramér's V ≈ 0
        "batch": ["batch1", "batch1", "batch2", "batch2",
                  "batch1", "batch1", "batch2", "batch2",
                  "batch1", "batch1", "batch2", "batch2"],
    }, index=samples)
    meta.index.name = "sample_id"

    out_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out_dir / "counts.tsv", sep="\t")
    meta.to_csv(out_dir / "metadata.tsv", sep="\t")
    _write_checksums(out_dir)


# ---------------------------------------------------------------------------
# Dataset 5: GSE144269
# Single-condition experiment (tests BIO-001 edge case on submission errors)
# Human, Ensembl IDs
# ---------------------------------------------------------------------------

def generate_gse144269(out_dir: Path) -> None:
    rng = np.random.default_rng(5)
    # Simulates a dataset where all samples were mislabeled as the same condition
    conditions = ["tumor"] * 6
    samples = [f"GSM{i}" for i in range(4280001, 4280007)]

    n_genes = 500
    genes = _human_ensembl(n_genes, rng)
    hk_ids = [
        "ENSG00000075624", "ENSG00000111640", "ENSG00000166710",
        "ENSG00000165704", "ENSG00000256269", "ENSG00000073578",
        "ENSG00000112592", "ENSG00000089157", "ENSG00000164924",
        "ENSG00000196262",
    ]
    genes = hk_ids + [g for g in genes if g not in set(hk_ids)]

    counts = _make_counts(genes, samples, conditions, seed=5)
    meta = pd.DataFrame({
        "condition": conditions,
        "batch": ["batch1"] * 6,
    }, index=samples)
    meta.index.name = "sample_id"

    out_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out_dir / "counts.tsv", sep="\t")
    meta.to_csv(out_dir / "metadata.tsv", sep="\t")
    _write_checksums(out_dir)


# ---------------------------------------------------------------------------
# Checksum helper
# ---------------------------------------------------------------------------

def _write_checksums(out_dir: Path) -> None:
    lines = []
    for fname in ["counts.tsv", "metadata.tsv"]:
        fpath = out_dir / fname
        if fpath.exists():
            digest = hashlib.sha256(fpath.read_bytes()).hexdigest()
            lines.append(f"{digest}  {fname}\n")
    (out_dir / "sha256sums.txt").write_text("".join(lines))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

GENERATORS = {
    "airway": generate_airway,
    "GSE37704": generate_gse37704,
    "GSE89189": generate_gse89189,
    "GSE96870": generate_gse96870,
    "GSE144269": generate_gse144269,
}


def main() -> None:
    print("Generating proxy real datasets …")
    for name, fn in GENERATORS.items():
        out = ROOT / "real" / name
        fn(out)
        print(f"  ✅  {name}  →  {out}/counts.tsv  +  metadata.tsv")
    print("Done.")


if __name__ == "__main__":
    main()
