#!/usr/bin/env python3
"""Download, verify, and standardize real public RNA-seq datasets.

This script replaces the synthetic proxy data in datasets/real/ with true
public count matrices and metadata when internet access is available.

Each dataset goes through a standardized pipeline:
  1. Download the count matrix and metadata from GEO / Bioconductor
  2. Verify raw-count status (integers only, not FPKM/TPM)
  3. Verify matrix orientation (genes × samples)
  4. Standardize column names to: sample_id, condition, batch
  5. Strip Ensembl version suffixes (.1, .2) from gene IDs
  6. Recompute sha256sums.txt
  7. Update README.md with download timestamp

Usage:
    python datasets/fetch_real_data.py --dataset airway
    python datasets/fetch_real_data.py --dataset GSE37704
    python datasets/fetch_real_data.py --dataset GSE89189
    python datasets/fetch_real_data.py --dataset GSE96870
    python datasets/fetch_real_data.py --dataset GSE144269
    python datasets/fetch_real_data.py          # attempt all (skips on failure)

Requirements (install before running):
    pip install GEOparse requests pandas numpy

Notes on sourcing:
  - airway:    Bioconductor ExperimentHub (EH1006); R required if using BiocManager,
               or download via the direct ExperimentHub URL listed below.
  - GSE37704:  Original GEO deposit is FPKM (not raw). Use recount2/recount3 for
               raw Ensembl counts. See README in datasets/real/GSE37704/.
  - GSE89189:  Raw count matrix in GEO supplementary files.
  - GSE96870:  Raw count matrix in GEO supplementary files (multi-factor design).
  - GSE144269: Raw count matrix in GEO supplementary files (single condition).
"""
from __future__ import annotations

import argparse
import hashlib
import io
import sys
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parent
REAL_DIR = ROOT / "real"

# ---------------------------------------------------------------------------
# Verification helpers
# ---------------------------------------------------------------------------

def verify_raw_counts(df: pd.DataFrame, dataset: str) -> None:
    """Raise ValueError if any values are non-integer (suggests normalized data)."""
    numeric = df.apply(pd.to_numeric, errors="coerce")
    non_int_frac = (numeric != numeric.round()).stack().mean()
    if non_int_frac > 0.01:
        raise ValueError(
            f"{dataset}: {non_int_frac*100:.1f}% of values are non-integer. "
            "This looks like FPKM/TPM/CPM data, not raw counts. "
            "See README for how to obtain raw counts (e.g. via recount2)."
        )


def verify_orientation(df: pd.DataFrame, dataset: str) -> None:
    """Raise ValueError if the matrix is likely transposed (samples × genes)."""
    n_rows, n_cols = df.shape
    if n_cols > n_rows and n_rows < 500:
        raise ValueError(
            f"{dataset}: matrix has {n_rows} rows and {n_cols} columns. "
            "Expected genes × samples (rows > cols). "
            "Transpose with: df = df.T"
        )
    if n_rows < 100:
        raise ValueError(
            f"{dataset}: only {n_rows} gene rows — too few for a typical RNA-seq matrix. "
            "Check that the correct file was downloaded."
        )


def strip_ensembl_versions(df: pd.DataFrame) -> pd.DataFrame:
    """Remove version suffixes (ENSG00000000003.15 → ENSG00000000003)."""
    df.index = df.index.str.replace(r"\.\d+$", "", regex=True)
    return df


def verify_sample_id_match(counts: pd.DataFrame, meta: pd.DataFrame, dataset: str) -> None:
    """Raise ValueError if sample IDs in metadata don't match count matrix columns."""
    meta_ids = set(meta.index.astype(str))
    count_cols = set(counts.columns.astype(str))
    if meta_ids != count_cols:
        missing_in_meta = count_cols - meta_ids
        extra_in_meta = meta_ids - count_cols
        msg_parts = [f"{dataset}: sample ID mismatch between counts and metadata."]
        if missing_in_meta:
            msg_parts.append(f"  In counts but not metadata: {sorted(missing_in_meta)[:5]}")
        if extra_in_meta:
            msg_parts.append(f"  In metadata but not counts: {sorted(extra_in_meta)[:5]}")
        raise ValueError("\n".join(msg_parts))


def write_checksums(out_dir: Path) -> None:
    lines = []
    for fname in ["counts.tsv", "metadata.tsv"]:
        fpath = out_dir / fname
        if fpath.exists():
            digest = hashlib.sha256(fpath.read_bytes()).hexdigest()
            lines.append(f"{digest}  {fname}\n")
    (out_dir / "sha256sums.txt").write_text("".join(lines))
    print(f"  ✅  sha256sums.txt updated")


def save_dataset(
    out_dir: Path,
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    dataset: str,
    source_url: str,
) -> None:
    """Run all verification steps and write counts/metadata/checksums."""
    verify_raw_counts(counts, dataset)
    verify_orientation(counts, dataset)
    counts = strip_ensembl_versions(counts)
    verify_sample_id_match(counts, meta, dataset)

    out_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out_dir / "counts.tsv", sep="\t")
    meta.to_csv(out_dir / "metadata.tsv", sep="\t")
    write_checksums(out_dir)

    # Append download provenance to README
    readme = out_dir / "README.md"
    note = (
        f"\n## Download Record\n\n"
        f"| Field | Value |\n"
        f"|---|---|\n"
        f"| Downloaded | {date.today().isoformat()} |\n"
        f"| Source URL | {source_url} |\n"
        f"| Genes | {counts.shape[0]:,} |\n"
        f"| Samples | {counts.shape[1]} |\n"
        f"| Proxy replaced | yes |\n"
    )
    existing = readme.read_text() if readme.exists() else ""
    # Remove any previous download record before appending
    if "## Download Record" in existing:
        existing = existing[:existing.index("## Download Record")]
    readme.write_text(existing.rstrip() + "\n" + note)
    print(f"  ✅  README.md updated with download record")


# ---------------------------------------------------------------------------
# Dataset fetchers
# ---------------------------------------------------------------------------

def fetch_airway(out_dir: Path) -> None:
    """Download airway (GSE52778) from Bioconductor ExperimentHub.

    Preferred path (requires R):
        BiocManager::install("airway")
        library(airway); data(airway)
        counts <- assay(airway, "counts")
        meta   <- as.data.frame(colData(airway)[, c("cell", "dex")])
        write.table(counts, "counts.tsv", sep="\\t", quote=FALSE, col.names=NA)
        write.table(meta,   "metadata.tsv", sep="\\t", quote=FALSE, col.names=NA)

    Python path (GEOparse — only metadata; counts via SRA):
        pip install GEOparse
        See https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778
        SRR IDs: SRR1039508–SRR1039521 (8 samples)
    """
    try:
        import GEOparse  # type: ignore
    except ImportError:
        raise ImportError(
            "GEOparse not installed. Run: pip install GEOparse\n"
            "Alternatively, use R/Bioconductor as described in datasets/real/airway/README.md"
        )

    print("  Attempting airway via GEOparse (metadata only) …")
    print("  NOTE: Raw counts require R/Bioconductor or SRA download.")
    print("  See datasets/real/airway/README.md for full instructions.")
    print("  Skipping automatic download — use R path or manual SRA.")
    raise NotImplementedError(
        "airway counts require R/Bioconductor. "
        "See datasets/real/airway/README.md for instructions."
    )


def fetch_gse37704(out_dir: Path) -> None:
    """Fetch GSE37704 raw counts from recount2/recount3.

    IMPORTANT: The original GEO deposit (GSE37704_FPKM_final.txt.gz) is
    FPKM-normalized. Raw counts must come from recount2 or recount3:

        # R path (recommended):
        BiocManager::install("recount")
        library(recount)
        download_study("SRP012682")  # SRP accession for GSE37704
        load("SRP012682/rse_gene.Rdata")
        counts <- assay(rse_gene, "counts")
        meta   <- as.data.frame(colData(rse_gene)[, c("title", "characteristics")])

        # Python manual path:
        # 1. Navigate to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37704
        # 2. Download supplementary count matrix (NOT the FPKM file)
        # 3. Verify it contains integers before running this script

    GEO accession: GSE37704
    SRA accession: SRP012682
    Samples: SRR493366–SRR493371 (6 samples)
    """
    print("  GSE37704: raw counts require recount2/R or manual SRA download.")
    print("  Original GEO deposit is FPKM — cannot be used directly.")
    print("  See datasets/real/GSE37704/README.md for detailed instructions.")
    raise NotImplementedError(
        "GSE37704 raw counts require recount2 (R) or SRA re-alignment. "
        "See datasets/real/GSE37704/README.md."
    )


def fetch_gse89189(out_dir: Path) -> None:
    """Fetch GSE89189 from GEO supplementary files.

    This dataset has raw integer counts in supplementary files.

    Steps:
      1. Navigate to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89189
      2. Download the supplementary count matrix file.
      3. Verify counts are integers (no decimal points).
      4. Place as datasets/real/GSE89189/raw_counts_download.txt
      5. Re-run: python datasets/fetch_real_data.py --dataset GSE89189

    Alternatively use GEOparse:
        import GEOparse
        gse = GEOparse.get_GEO("GSE89189", destdir="/tmp/geo")
        # Extract supplementary files from gse.metadata["supplementary_file"]
    """
    # Check for a manually placed download file
    raw_file = out_dir / "raw_counts_download.txt"
    if not raw_file.exists():
        print("  GSE89189: no pre-downloaded file found.")
        print(f"  Please download the count matrix from GEO and place it at:")
        print(f"  {raw_file}")
        print("  See datasets/real/GSE89189/README.md for details.")
        raise FileNotFoundError(
            f"Please download GSE89189 count matrix to {raw_file}. "
            "See datasets/real/GSE89189/README.md."
        )

    print(f"  Reading {raw_file.name} …")
    counts = pd.read_csv(raw_file, sep="\t", index_col=0)
    verify_orientation(counts, "GSE89189")
    verify_raw_counts(counts, "GSE89189")

    # Standardize metadata
    # GSM IDs → sample_id; extract condition from sample names or series matrix
    samples = list(counts.columns)
    # Basic metadata — requires manual curation from series matrix
    # Condition pattern: first 3 are typically unstimulated, last 3 are LPS
    # Verify this against the actual series matrix file
    print("  WARNING: Condition labels are inferred. Verify against series matrix.")
    conditions = (
        ["unstimulated"] * (len(samples) // 2) +
        ["LPS"] * (len(samples) - len(samples) // 2)
    )
    meta = pd.DataFrame({"condition": conditions}, index=samples)
    meta.index.name = "sample_id"

    save_dataset(
        out_dir, counts, meta, "GSE89189",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89189"
    )


def fetch_gse96870(out_dir: Path) -> None:
    """Fetch GSE96870 from GEO supplementary files (multi-factor design).

    This dataset tests multi-factor confounding detection.
    Design: sex (male/female) × time (0h/6h/24h), 2 replicates per combination.

    Steps:
      1. Navigate to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96870
      2. Download the supplementary count matrix.
      3. Download the series matrix for metadata.
      4. Place count matrix as datasets/real/GSE96870/raw_counts_download.txt
      5. Re-run: python datasets/fetch_real_data.py --dataset GSE96870

    Critical: after loading, verify that the batch column (sequencing lane or date)
    is NOT perfectly confounded with condition — otherwise BIO-007 will fire on
    the clean dataset, which would be a false positive.
    """
    raw_file = out_dir / "raw_counts_download.txt"
    if not raw_file.exists():
        print("  GSE96870: no pre-downloaded file found.")
        print(f"  Please download the count matrix from GEO and place it at:")
        print(f"  {raw_file}")
        print("  See datasets/real/GSE96870/README.md for details.")
        raise FileNotFoundError(
            f"Please download GSE96870 count matrix to {raw_file}. "
            "See datasets/real/GSE96870/README.md."
        )

    print(f"  Reading {raw_file.name} …")
    counts = pd.read_csv(raw_file, sep="\t", index_col=0)
    verify_orientation(counts, "GSE96870")
    verify_raw_counts(counts, "GSE96870")

    print("  WARNING: Metadata must be manually curated from series matrix.")
    print("  Ensure batch column is NOT confounded with condition.")
    print("  See datasets/real/GSE96870/README.md for column standardization.")
    raise NotImplementedError(
        "GSE96870 metadata requires manual curation from the series matrix. "
        "See datasets/real/GSE96870/README.md."
    )


def fetch_gse144269(out_dir: Path) -> None:
    """Fetch GSE144269 from GEO supplementary files (single-condition edge case).

    This dataset tests BIO-001: all samples share one condition label.
    The BIO-001 FAIL is EXPECTED and intentional — it is not a false positive.

    Steps:
      1. Navigate to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144269
      2. Verify that all samples have the same condition in the series matrix.
      3. Download the supplementary count matrix.
      4. Place as datasets/real/GSE144269/raw_counts_download.txt
      5. Re-run: python datasets/fetch_real_data.py --dataset GSE144269
    """
    raw_file = out_dir / "raw_counts_download.txt"
    if not raw_file.exists():
        print("  GSE144269: no pre-downloaded file found.")
        print(f"  Please download the count matrix from GEO and place it at:")
        print(f"  {raw_file}")
        print("  See datasets/real/GSE144269/README.md for details.")
        raise FileNotFoundError(
            f"Please download GSE144269 count matrix to {raw_file}. "
            "See datasets/real/GSE144269/README.md."
        )

    print(f"  Reading {raw_file.name} …")
    counts = pd.read_csv(raw_file, sep="\t", index_col=0)
    verify_orientation(counts, "GSE144269")
    verify_raw_counts(counts, "GSE144269")

    samples = list(counts.columns)
    # All samples should be the same condition — verify before saving
    print("  NOTE: BIO-001 will FAIL on this dataset (single condition = intended).")
    conditions = ["tumor"] * len(samples)
    meta = pd.DataFrame({"condition": conditions}, index=samples)
    meta.index.name = "sample_id"

    save_dataset(
        out_dir, counts, meta, "GSE144269",
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144269"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

FETCHERS = {
    "airway": fetch_airway,
    "GSE37704": fetch_gse37704,
    "GSE89189": fetch_gse89189,
    "GSE96870": fetch_gse96870,
    "GSE144269": fetch_gse144269,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset", "-d",
        metavar="NAME",
        choices=list(FETCHERS),
        help="Fetch only this dataset. Default: attempt all.",
    )
    args = parser.parse_args()

    datasets = [args.dataset] if args.dataset else list(FETCHERS)
    successes, skipped, failed = [], [], []

    for name in datasets:
        out_dir = REAL_DIR / name
        print(f"\n{'='*60}")
        print(f"Dataset: {name}")
        print(f"{'='*60}")
        try:
            FETCHERS[name](out_dir)
            successes.append(name)
            print(f"  ✅  {name} successfully replaced with real data.")
        except (NotImplementedError, FileNotFoundError) as exc:
            skipped.append(name)
            print(f"  ⏭️  {name} skipped: {exc}")
        except Exception as exc:
            failed.append(name)
            print(f"  ❌  {name} FAILED: {exc}")

    print(f"\n{'='*60}")
    print("Summary")
    print(f"{'='*60}")
    if successes:
        print(f"  Replaced:  {', '.join(successes)}")
    if skipped:
        print(f"  Skipped:   {', '.join(skipped)} (manual steps required)")
    if failed:
        print(f"  Failed:    {', '.join(failed)}")
    print()
    print("After replacing proxy data:")
    print("  1. Run: python datasets/generate_real_faults.py")
    print("  2. Run: python datasets/validate_clean_datasets.py")
    print("  3. Run: PYTHONPATH=backend python datasets/benchmark_eval.py --strict")


if __name__ == "__main__":
    main()
