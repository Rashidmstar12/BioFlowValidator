#!/usr/bin/env python3
"""Manual real-data ingestion workflow for BioFlowValidator.

This script accepts locally placed count matrices and metadata files, verifies
their structure, standardizes column names, computes SHA-256 checksums, and
writes provenance records — enabling true-public-data benchmarking outside CI.

Usage (single dataset):
    python datasets/ingest_real_data.py \\
        --dataset GSE60450 \\
        --counts /path/to/GSE60450_Lactation-GenewiseCounts.txt \\
        --metadata /path/to/GSE60450_metadata.csv

Usage (verify an already-placed dataset, no new files):
    python datasets/ingest_real_data.py --dataset GSE60450 --verify-only

Workflow
--------
1.  Read the input count matrix and metadata files.
2.  Verify count matrix:
    - integer values (NRM-001 pre-check)
    - genes × samples orientation (rows >> columns)
    - no duplicate gene IDs
3.  Standardize metadata:
    - rename columns to sample_id, condition, batch (best-effort)
    - ensure sample_ids match count matrix column headers
4.  Write standardized files to datasets/real/<DATASET>/:
    - counts.tsv
    - metadata.tsv
    - provenance.json  (source file, download date, checksums, operator notes)
    - sha256sums.txt
5.  Update benchmark_manifest.yaml with checksums and source_type=true_public_data.

This script intentionally does NOT download any files.  Place downloaded files
locally first.  See datasets/real/<DATASET>/README.md for download instructions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).parent
REAL_DIR = ROOT / "real"
BACKEND = ROOT.parent / "backend"
sys.path.insert(0, str(BACKEND))

# Known column name aliases for standardization
_CONDITION_ALIASES = {
    "treatment", "group", "condition", "phenotype", "cell_type", "celltype",
    "cell.type", "tissue", "genotype", "stimulation",
}
_BATCH_ALIASES = {
    "batch", "lane", "plate", "run", "flow_cell", "flowcell",
    "sequencing_batch", "replicate", "donor", "patient", "subject",
    "individual", "sample_batch",
}
_SAMPLE_ID_ALIASES = {
    "sample_id", "sampleid", "sample", "sample.id", "gsm", "srr",
    "run_accession", "biosample",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_checksums(out_dir: Path) -> dict[str, str]:
    sums: dict[str, str] = {}
    lines = []
    for fname in ["counts.tsv", "metadata.tsv"]:
        fpath = out_dir / fname
        if fpath.exists():
            digest = _sha256(fpath)
            sums[fname] = digest
            lines.append(f"{digest}  {fname}\n")
    (out_dir / "sha256sums.txt").write_text("".join(lines))
    return sums


def _read_counts(path: Path) -> pd.DataFrame:
    """Read a count matrix from TSV, CSV, or whitespace-delimited file.

    Special handling: if the first non-index column is named exactly ``Length``
    (case-sensitive), it is treated as a gene-annotation column and dropped.
    This matches the layout of GSE60450_Lactation-GenewiseCounts.txt, where GEO
    places a gene-length annotation as the second column before the sample counts.
    """
    suffix = path.suffix.lower()
    sep: str
    if suffix in {".tsv", ".txt"}:
        sep = "\t"
    elif suffix == ".csv":
        sep = ","
    else:
        # Try to detect
        first_line = path.read_text(errors="replace").split("\n")[0]
        sep = "\t" if "\t" in first_line else ","

    df = pd.read_csv(path, sep=sep, index_col=0)

    # Drop gene-annotation column "Length" if present (e.g. GSE60450 raw GEO file).
    if "Length" in df.columns:
        print(
            "  ℹ  Dropping 'Length' annotation column "
            f"(found in {path.name}; not a sample column)."
        )
        df = df.drop(columns=["Length"])

    return df


def _standardize_metadata(meta: pd.DataFrame, dataset: str) -> tuple[pd.DataFrame, list[str]]:
    """Rename columns to sample_id, condition, batch where possible.

    Returns (standardized_df, warnings).
    """
    warnings: list[str] = []
    col_map: dict[str, str] = {}

    # Normalize column names to lowercase for matching
    for col in meta.columns:
        col_lower = col.lower().replace("-", "_").replace(" ", "_")
        if col_lower in _CONDITION_ALIASES and "condition" not in col_map.values():
            col_map[col] = "condition"
        elif col_lower in _BATCH_ALIASES and "batch" not in col_map.values():
            col_map[col] = "batch"

    if col_map:
        meta = meta.rename(columns=col_map)

    # Ensure index is named sample_id
    if meta.index.name is None or meta.index.name.lower() not in {"sample_id", "sample", "gsm", "srr"}:
        # Check if a column looks like sample IDs
        for col in meta.columns:
            col_lower = col.lower().replace("-", "_").replace(" ", "_")
            if col_lower in _SAMPLE_ID_ALIASES:
                meta = meta.set_index(col)
                break
    meta.index.name = "sample_id"

    if "condition" not in meta.columns:
        warnings.append(
            "No 'condition' column found or inferred. "
            "Please add a 'condition' column manually before running the benchmark."
        )
    if "batch" not in meta.columns:
        warnings.append(
            "No 'batch' column found or inferred. "
            "BIO-007 (batch confounding) will SKIP without a batch column."
        )
    return meta, warnings


def _verify_counts(df: pd.DataFrame) -> list[str]:
    """Return a list of issues found in the count matrix.  Empty list = clean."""
    issues = []
    n_rows, n_cols = df.shape

    # Orientation check
    if n_cols > n_rows:
        issues.append(
            f"Matrix appears TRANSPOSED: {n_rows} rows × {n_cols} columns "
            "(expected genes × samples, i.e. rows >> columns). "
            "Transpose before ingestion."
        )

    # Integer check
    numeric = df.apply(pd.to_numeric, errors="coerce")
    non_int = int(numeric.notna().sum().sum()) - int((numeric % 1 == 0).sum().sum())
    if non_int > 0:
        frac = non_int / max(numeric.notna().sum().sum(), 1)
        issues.append(
            f"{non_int} values ({frac*100:.1f}%) are non-integer. "
            "Use the raw-count matrix, not TPM/CPM/RPKM."
        )

    # Duplicate gene IDs
    dup = df.index.duplicated().sum()
    if dup > 0:
        issues.append(f"{dup} duplicate gene IDs found. Remove or deduplicate before ingestion.")

    # Negative values
    if (numeric < 0).any().any():
        issues.append("Negative count values found. This is not a valid raw-count matrix.")

    return issues


def ingest(
    dataset: str,
    counts_path: Path | None,
    meta_path: Path | None,
    notes: str = "",
    verify_only: bool = False,
) -> bool:
    """Ingest and standardize a real dataset.

    Returns True on success, False on blocking error.
    """
    out_dir = REAL_DIR / dataset
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  Dataset: {dataset}")
    print(f"{'='*60}")

    if verify_only:
        counts_path = out_dir / "counts.tsv"
        meta_path = out_dir / "metadata.tsv"
        if not counts_path.exists():
            print(f"  ❌ {counts_path} not found. Use --counts to provide the file.")
            return False
        print("  Mode: verify-only (using existing files)")

    # --- Read counts ---
    if counts_path is None or not counts_path.exists():
        print(f"  ❌ Counts file not found: {counts_path}")
        return False

    print(f"  Reading counts: {counts_path}")
    try:
        counts = _read_counts(counts_path)
    except Exception as exc:
        print(f"  ❌ Failed to read counts: {exc}")
        return False

    print(f"  Shape: {counts.shape[0]:,} genes × {counts.shape[1]:,} samples")

    # --- Verify counts ---
    issues = _verify_counts(counts)
    if issues:
        print("  ⚠  Count matrix issues:")
        for issue in issues:
            print(f"     • {issue}")
        blocking = [i for i in issues if "TRANSPOSED" in i or "non-integer" in i or "Negative" in i]
        if blocking:
            print("  ❌ Blocking issues found — resolve before ingestion.")
            return False
    else:
        print("  ✅ Count matrix verified: integer values, correct orientation, no duplicates.")

    # --- Read metadata ---
    meta: pd.DataFrame | None = None
    meta_warnings: list[str] = []
    if meta_path is not None and meta_path.exists():
        print(f"  Reading metadata: {meta_path}")
        try:
            suffix = meta_path.suffix.lower()
            sep = "\t" if suffix in {".tsv", ".txt"} else ","
            meta = pd.read_csv(meta_path, sep=sep, index_col=0)
            meta, meta_warnings = _standardize_metadata(meta, dataset)
            # Verify sample ID alignment
            count_cols = set(counts.columns.astype(str))
            meta_idx = set(meta.index.astype(str))
            overlap = count_cols & meta_idx
            only_counts = count_cols - meta_idx
            only_meta = meta_idx - count_cols
            print(
                f"  Metadata: {len(meta)} samples; "
                f"{len(overlap)} matched, {len(only_counts)} only in counts, "
                f"{len(only_meta)} only in metadata."
            )
            if only_counts:
                meta_warnings.append(
                    f"{len(only_counts)} count samples not in metadata: "
                    f"{sorted(only_counts)[:5]}…"
                )
        except Exception as exc:
            print(f"  ⚠  Failed to read metadata: {exc}. Proceeding without metadata.")
            meta = None
    else:
        meta_warnings.append("No metadata file provided. BIO-007, SMP-004, and BIO-001 will SKIP.")
        print("  ℹ  No metadata provided.")

    if meta_warnings:
        print("  ⚠  Metadata warnings:")
        for w in meta_warnings:
            print(f"     • {w}")

    if verify_only:
        print("  ✅ Verification complete.")
        return True

    # --- Write standardized files ---
    print(f"  Writing to: {out_dir}/")
    counts.to_csv(out_dir / "counts.tsv", sep="\t")
    if meta is not None:
        meta.to_csv(out_dir / "metadata.tsv", sep="\t")

    sums = _write_checksums(out_dir)

    # Mark this dataset as true_public_data
    (out_dir / "benchmark_mode.txt").write_text("true_public_data\n")

    # --- Write provenance record ---
    provenance: dict[str, Any] = {
        "dataset": dataset,
        "source_type": "true_public_data",
        "ingested_at": datetime.now(timezone.utc).isoformat(),
        "source_counts_file": str(counts_path),
        "source_metadata_file": str(meta_path) if meta_path else None,
        "shape": {"n_genes": counts.shape[0], "n_samples": counts.shape[1]},
        "checksums": sums,
        "count_matrix_issues": issues,
        "metadata_warnings": meta_warnings,
        "operator_notes": notes,
    }
    (out_dir / "provenance.json").write_text(json.dumps(provenance, indent=2))

    print(f"  ✅ Ingested: {counts.shape[0]:,} genes × {counts.shape[1]:,} samples")
    print(f"     SHA-256 counts:   {sums.get('counts.tsv', 'N/A')[:16]}…")
    print(f"     SHA-256 metadata: {sums.get('metadata.tsv', 'N/A')[:16]}…")
    print(f"     Provenance:       {out_dir}/provenance.json")
    print()
    print("  Next steps:")
    print(f"    1. Review {out_dir}/metadata.tsv — check condition and batch columns.")
    print(f"    2. Update datasets/benchmark_manifest.yaml:")
    print(f"       source_type: true_public_data")
    print(f"       counts_sha256: {sums.get('counts.tsv', '')}")
    print(f"       metadata_sha256: {sums.get('metadata.tsv', '')}")
    print("    3. Run: PYTHONPATH=backend python datasets/validate_clean_datasets.py")
    print("    4. Run: python datasets/generate_real_faults.py --dataset", dataset)
    print("    5. Run: PYTHONPATH=backend python datasets/benchmark_eval.py --strict")

    return True


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--dataset", required=True, help="Dataset name (e.g. GSE60450).")
    parser.add_argument("--counts", type=Path, default=None,
                        help="Path to the raw count matrix file (TSV, CSV, or TXT).")
    parser.add_argument("--metadata", type=Path, default=None,
                        help="Path to the sample metadata file (TSV or CSV).")
    parser.add_argument("--notes", default="",
                        help="Free-text notes to record in provenance.json.")
    parser.add_argument("--verify-only", action="store_true",
                        help="Verify existing files in datasets/real/<DATASET>/ without re-ingesting.")
    args = parser.parse_args()

    success = ingest(
        dataset=args.dataset,
        counts_path=args.counts,
        meta_path=args.metadata,
        notes=args.notes,
        verify_only=args.verify_only,
    )
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
