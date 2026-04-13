#!/usr/bin/env python3
"""Fault injection script for Phase 8 benchmark datasets.

For each clean real (or proxy) dataset in datasets/real/, this script generates
7 faulty variants and writes them to datasets/real_faulty/<dataset>/.
Each variant introduces exactly one injected fault, along with a ground-truth
YAML file describing the expected validator output.

Fault types:
  F1  sample_label_mismatch   — rename 2 sample IDs in metadata
  F2  mixed_gene_ids          — replace first 50 gene IDs with HGNC symbols
  F3  transposed_matrix       — swap rows and columns of count matrix
  F4  batch_confounding       — make batch perfectly confounded with condition
  F5  low_replicates          — drop all but 1 sample from one condition
  F6  library_size_imbalance  — multiply one sample's counts by 50×
  F7  normalized_not_raw      — convert counts to CPM (floating-point values)

Usage:
    python datasets/generate_real_faults.py [--dataset NAME]
    python datasets/generate_real_faults.py           # process all datasets
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).parent
REAL_DIR = ROOT / "real"
FAULTY_DIR = ROOT / "real_faulty"

# ---------------------------------------------------------------------------
# Representative HGNC symbols used for F2 (mixed gene IDs)
# These are well-known human gene symbols that the _classify_id function
# will correctly classify as "symbol" type.
# ---------------------------------------------------------------------------
_HGNC_SYMBOLS = [
    "TP53", "ACTB", "GAPDH", "MYC", "BRCA1", "EGFR", "KRAS", "PTEN",
    "RB1", "APC", "CDH1", "VHL", "SMAD4", "PIK3CA", "BRAF", "ALK",
    "FGFR1", "NRAS", "IDH1", "CDKN2A", "MDM2", "BCL2", "VEGFA",
    "FLT3", "NPM1", "DNMT3A", "TET2", "RUNX1", "CEBPA", "NF1",
    "PTPN11", "JAK2", "STAT3", "KIT", "PDGFRA", "ABL1", "BCR",
    "PML", "RARA", "NOTCH1", "FBXW7", "SF3B1", "SRSF2", "U2AF1",
    "ASXL1", "EZH2", "SUZ12", "KDM6A", "ARID1A", "CREBBP",
]


def _load_dataset(dataset_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load counts.tsv and metadata.tsv from a dataset directory."""
    counts = pd.read_csv(dataset_dir / "counts.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(dataset_dir / "metadata.tsv", sep="\t", index_col=0)
    return counts, meta


def _save_variant(
    out_dir: Path,
    fault_type: str,
    counts: pd.DataFrame,
    meta: pd.DataFrame | None,
    truth: dict,
) -> None:
    """Write counts, metadata, and ground-truth YAML for one fault variant."""
    out_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out_dir / f"{fault_type}_counts.tsv", sep="\t")
    if meta is not None:
        meta.to_csv(out_dir / f"{fault_type}_meta.tsv", sep="\t")
    (out_dir / f"{fault_type}.yaml").write_text(
        yaml.dump(truth, default_flow_style=False, sort_keys=False)
    )


def _find_condition_col(meta: pd.DataFrame) -> str | None:
    _CONDITION_COLS = {"condition", "group", "treatment", "genotype", "cell_type", "dex"}
    for col in meta.columns:
        if col.lower() in _CONDITION_COLS:
            return col
    for col in meta.columns:
        if meta[col].dtype == object:
            return col
    return None


# ---------------------------------------------------------------------------
# Fault injectors
# ---------------------------------------------------------------------------

def inject_sample_label_mismatch(
    counts: pd.DataFrame, meta: pd.DataFrame, dataset: str
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """F1: Rename 2 sample IDs in metadata so they no longer match the count matrix."""
    new_meta = meta.copy()
    samples_to_rename = list(meta.index[:2])
    new_index = list(meta.index)
    for i, s in enumerate(samples_to_rename):
        new_index[i] = f"{s}_WRONG"
    new_meta.index = new_index

    truth = {
        "dataset": dataset,
        "fault_type": "sample_label_mismatch",
        "description": (
            f"2 sample IDs in metadata renamed with _WRONG suffix: {samples_to_rename}"
        ),
        "expected_fails": [
            {"rule_id": "SMP-001", "severity": "ERROR", "min_affected_items": 2}
        ],
        "expected_passes": [],
        "expected_skips": [],
        "notes": (
            "SMP-001 must FAIL because 2 metadata samples are not in the count matrix "
            "and 2 count matrix samples are not in the metadata."
        ),
    }
    return counts, new_meta, truth


def inject_mixed_gene_ids(
    counts: pd.DataFrame, meta: pd.DataFrame, dataset: str
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """F2: Replace first 50 non-housekeeping gene IDs with HGNC gene symbols."""
    from app.rules.biology.rules import _HUMAN_HK_ENSEMBL, _MOUSE_HK_ENSEMBL

    # Collect all housekeeping IDs (human + mouse) to skip over them
    hk_ensembl = set(_HUMAN_HK_ENSEMBL.keys()) | set(_MOUSE_HK_ENSEMBL.keys())

    new_index = list(counts.index)
    n_replaced = 0
    symbols = iter(_HGNC_SYMBOLS)
    for pos, gene_id in enumerate(new_index):
        if n_replaced >= 50:
            break
        # Skip housekeeping gene rows so BIO-006 doesn't fire as a cascade
        if str(gene_id).split(".")[0].upper() in {h.upper() for h in hk_ensembl}:
            continue
        try:
            new_index[pos] = next(symbols)
            n_replaced += 1
        except StopIteration:
            break

    new_counts = counts.copy()
    new_counts.index = new_index

    truth = {
        "dataset": dataset,
        "fault_type": "mixed_gene_ids",
        "description": f"{n_replaced} non-housekeeping gene IDs replaced with HGNC gene symbols",
        "expected_fails": [
            {"rule_id": "GEN-001", "severity": "ERROR", "min_affected_items": 1}
        ],
        "expected_passes": [],
        "expected_skips": [],
        "notes": (
            "GEN-001 must FAIL because HGNC symbols and Ensembl IDs are mixed. "
            "Housekeeping gene rows are preserved so BIO-006 does not fire as a cascade."
        ),
    }
    return new_counts, meta, truth


def inject_transposed_matrix(
    counts: pd.DataFrame, meta: pd.DataFrame, dataset: str
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """F3: Transpose the count matrix (samples × genes instead of genes × samples)."""
    transposed = counts.T

    truth = {
        "dataset": dataset,
        "fault_type": "transposed_matrix",
        "description": "Count matrix transposed: now samples × genes instead of genes × samples",
        "expected_fails": [
            {"rule_id": "FMT-008", "severity": "WARNING", "min_affected_items": 1}
        ],
        "expected_passes": [],
        "expected_skips": [],
        "notes": (
            "FMT-008 must FAIL because n_rows < n_cols and n_rows < 500. "
            "Other rules may also fire because sample names look like gene IDs."
        ),
    }
    return transposed, meta, truth


def inject_batch_confounding(
    counts: pd.DataFrame, meta: pd.DataFrame, dataset: str
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """F4: Make batch perfectly confounded with condition."""
    new_meta = meta.copy()
    cond_col = _find_condition_col(new_meta)

    if cond_col is None:
        # Should not happen for proxy datasets
        raise ValueError(f"No condition column found in {dataset} metadata")

    conditions = new_meta[cond_col].values
    unique_conds = list(dict.fromkeys(conditions))  # preserve order

    if len(unique_conds) < 2:
        # Single-condition dataset (e.g. GSE144269): confounding test requires ≥2 conditions.
        # Add a batch column anyway, but BIO-007 will SKIP because V cannot be computed.
        new_meta["batch"] = [f"batch_{chr(65 + (i % 2))}" for i in range(len(conditions))]
        expected_fails: list[dict] = []
        primary_note = (
            f"Dataset '{dataset}' has only one condition; batch-condition confounding "
            "cannot be computed (Cramér's V requires ≥2 categories). "
            "BIO-007 will SKIP or PASS. BIO-001 remains the primary expected fail."
        )
    else:
        # Assign each condition to its own batch → perfect confounding
        cond_to_batch = {c: f"batch_{chr(65+i)}" for i, c in enumerate(unique_conds)}
        new_meta["batch"] = [cond_to_batch[c] for c in conditions]
        expected_fails = [
            {"rule_id": "BIO-007", "severity": "ERROR", "min_affected_items": 1}
        ]
        primary_note = (
            f"Batch column added so each condition maps to exactly one batch "
            f"(Cramér's V = 1.0). Condition column: '{cond_col}'"
        )

    truth = {
        "dataset": dataset,
        "fault_type": "batch_confounding",
        "description": primary_note,
        "expected_fails": expected_fails,
        "expected_passes": [],
        "expected_skips": [],
        "notes": primary_note,
    }
    return counts, new_meta, truth


def inject_low_replicates(
    counts: pd.DataFrame, meta: pd.DataFrame, dataset: str
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """F5: Drop all but 1 sample from one condition group."""
    cond_col = _find_condition_col(meta)
    if cond_col is None:
        raise ValueError(f"No condition column found in {dataset} metadata")

    conditions = meta[cond_col].value_counts()
    # Pick the condition with the most replicates to drop from
    target_cond = conditions.idxmax()
    target_samples = meta[meta[cond_col] == target_cond].index.tolist()
    # Keep only the first sample from this condition
    samples_to_drop = target_samples[1:]

    new_meta = meta.drop(index=samples_to_drop)
    new_counts = counts.drop(columns=samples_to_drop)

    truth = {
        "dataset": dataset,
        "fault_type": "low_replicates",
        "description": (
            f"Dropped {len(samples_to_drop)} sample(s) from condition '{target_cond}', "
            f"leaving only 1 replicate"
        ),
        "expected_fails": [
            {"rule_id": "SMP-004", "severity": "ERROR", "min_affected_items": 1}
        ],
        "expected_passes": [],
        "expected_skips": [],
        "notes": (
            f"SMP-004 must FAIL/ERROR because condition '{target_cond}' "
            "has only 1 replicate (< 2 required)."
        ),
    }
    return new_counts, new_meta, truth


def inject_library_size_imbalance(
    counts: pd.DataFrame, meta: pd.DataFrame, dataset: str
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """F6: Multiply one sample's counts by 50× to create extreme library size imbalance."""
    new_counts = counts.copy()
    target_sample = counts.columns[0]
    new_counts[target_sample] = (counts[target_sample] * 50).astype(int)

    original_size = int(counts[target_sample].sum())
    new_size = int(new_counts[target_sample].sum())

    truth = {
        "dataset": dataset,
        "fault_type": "library_size_imbalance",
        "description": (
            f"Sample '{target_sample}' counts multiplied by 50× "
            f"(original: {original_size:,}, new: {new_size:,})"
        ),
        "expected_fails": [
            {"rule_id": "NRM-002", "severity": "WARNING", "min_affected_items": 1}
        ],
        "expected_passes": [],
        "expected_skips": [],
        "notes": (
            "NRM-002 must FAIL because max/min library size ratio will exceed 10×. "
            f"The inflated sample is '{target_sample}'."
        ),
    }
    return new_counts, meta, truth


def inject_normalized_not_raw(
    counts: pd.DataFrame, meta: pd.DataFrame, dataset: str
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """F7: Convert all counts to CPM (non-integer floating-point values)."""
    lib_sizes = counts.sum(axis=0)
    # CPM = count / lib_size * 1e6
    cpm = counts.div(lib_sizes, axis=1) * 1_000_000

    truth = {
        "dataset": dataset,
        "fault_type": "normalized_not_raw",
        "description": "All count matrix values converted to CPM (floating-point, non-integer)",
        "expected_fails": [
            {"rule_id": "NRM-001", "severity": "ERROR", "min_affected_items": 0}
        ],
        "expected_passes": [],
        "expected_skips": [],
        "notes": (
            "NRM-001 must FAIL/ERROR because > 1% of values are non-integer. "
            "CPM values are floating-point; the fraction non-integer will be ~1.0."
        ),
    }
    return cpm, meta, truth


# ---------------------------------------------------------------------------
# Dataset-level baseline expected failures
# These rules are expected to fire on ALL variants of the given dataset,
# regardless of which fault is injected.  They are merged into every
# fault variant's ground-truth YAML.
# ---------------------------------------------------------------------------

_DATASET_BASELINE_FAILS: dict[str, list[dict]] = {
    # GSE144269: all samples labeled as the same condition, so BIO-001 always fires
    "GSE144269": [
        {"rule_id": "BIO-001", "severity": "ERROR", "min_affected_items": 0}
    ],
}

# ---------------------------------------------------------------------------
# Fault type registry
# ---------------------------------------------------------------------------

FAULT_INJECTORS = {
    "sample_label_mismatch": inject_sample_label_mismatch,
    "mixed_gene_ids": inject_mixed_gene_ids,
    "transposed_matrix": inject_transposed_matrix,
    "batch_confounding": inject_batch_confounding,
    "low_replicates": inject_low_replicates,
    "library_size_imbalance": inject_library_size_imbalance,
    "normalized_not_raw": inject_normalized_not_raw,
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def generate_faults_for_dataset(dataset_name: str) -> None:
    src_dir = REAL_DIR / dataset_name
    if not src_dir.exists():
        print(f"  ⚠️  {dataset_name}: source directory not found, skipping")
        return

    counts_path = src_dir / "counts.tsv"
    meta_path = src_dir / "metadata.tsv"
    if not counts_path.exists():
        print(f"  ⚠️  {dataset_name}: counts.tsv not found, skipping")
        return

    counts, meta = _load_dataset(src_dir)
    out_dir = FAULTY_DIR / dataset_name
    baseline_fails = _DATASET_BASELINE_FAILS.get(dataset_name, [])
    print(f"  Processing {dataset_name} ({counts.shape[0]} genes × {counts.shape[1]} samples)…")

    for fault_type, injector in FAULT_INJECTORS.items():
        try:
            faulty_counts, faulty_meta, truth = injector(counts.copy(), meta.copy(), dataset_name)
            # Merge baseline expected fails (preserve injector's primary fail first)
            primary_ids = {f["rule_id"] for f in truth["expected_fails"]}
            for b in baseline_fails:
                if b["rule_id"] not in primary_ids:
                    truth["expected_fails"].append(b)
            _save_variant(out_dir, fault_type, faulty_counts, faulty_meta, truth)
            print(f"    ✅  {fault_type}")
        except Exception as exc:
            print(f"    ❌  {fault_type}: {exc}")


def main(datasets: list[str] | None = None) -> None:
    if datasets is None:
        datasets = [d.name for d in REAL_DIR.iterdir() if d.is_dir()]
        datasets.sort()

    print(f"Generating fault variants for {len(datasets)} dataset(s) …")
    for ds in datasets:
        generate_faults_for_dataset(ds)
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset", "-d",
        metavar="NAME",
        help="Process only this dataset (e.g. airway). Default: all.",
    )
    args = parser.parse_args()
    main([args.dataset] if args.dataset else None)
