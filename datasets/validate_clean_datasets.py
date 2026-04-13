#!/usr/bin/env python3
"""Smoke test: run the validator on every clean real dataset and classify results.

This script implements the Phase A / Phase B validation from the real-data
replacement plan:

  - Runs the validator on each datasets/real/<dataset>/ directory
  - Classifies every non-PASS result as one of:
      EXPECTED   — known intentional fail (e.g. BIO-001 on GSE144269)
      WARNING    — rule fires but is a WARNING-severity (may be acceptable)
      FALSE_POS  — ERROR-severity rule fires on a clean dataset
  - Exits 0 if 0 FALSE_POS found, exits 1 otherwise (suitable for CI gate)

Usage:
    PYTHONPATH=backend python datasets/validate_clean_datasets.py
    PYTHONPATH=backend python datasets/validate_clean_datasets.py --verbose
    PYTHONPATH=backend python datasets/validate_clean_datasets.py --strict

Options:
    --verbose   Print all PASS and SKIP results in addition to FAIL/ERROR
    --strict    Treat WARNING-severity FP as failures (default: only ERROR)
    --dataset   Run only one dataset (e.g. --dataset airway)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).parent
BACKEND = ROOT.parent / "backend"
sys.path.insert(0, str(BACKEND))

from app.engine.parser import parse_files
from app.engine.runner import run_all
from app.models.rule_result import RuleResult

REAL_DIR = ROOT / "real"

# ---------------------------------------------------------------------------
# Known intentional failures on clean real datasets.
# These are NOT false positives — they represent the intended scientific test
# for each dataset.  Add any dataset-specific expected rules here.
# ---------------------------------------------------------------------------
EXPECTED_FAILS: dict[str, list[str]] = {
    # GSE144269: all samples share one condition → BIO-001 is expected
    "GSE144269": ["BIO-001"],
}

# ---------------------------------------------------------------------------
# Threshold documentation — for paper Table S1
# ---------------------------------------------------------------------------
RULE_THRESHOLDS = {
    "NRM-002": "Library size max/min > 10× (Conesa et al. 2016)",
    "BIO-005": "MT fraction > 30% per sample (conservative for bulk; Ilicic et al. 2016)",
    "BIO-007": "Cramér's V ≥ 0.999 (perfect confounding); V > 0.7 (near-perfect)",
    "SMP-005": "log1p-Pearson r ≥ 0.999 between sample pairs",
    "GEN-001": "Mixed namespaces detected across first 2000 gene IDs",
    "FMT-008": "n_cols > n_rows AND n_rows < 500 (transposed matrix)",
}


def run_clean_dataset(
    dataset: str,
    verbose: bool = False,
    strict: bool = False,
) -> tuple[int, int, int]:
    """Run validator on a clean dataset.

    Returns
    -------
    (n_expected, n_warning_fp, n_error_fp)
    """
    ds_dir = REAL_DIR / dataset
    counts_path = ds_dir / "counts.tsv"
    meta_path = ds_dir / "metadata.tsv"

    if not counts_path.exists():
        print(f"  {dataset}: counts.tsv not found — skipping")
        return 0, 0, 0

    count_bytes = counts_path.read_bytes()
    meta_bytes = meta_path.read_bytes() if meta_path.exists() else None
    ctx = parse_files(
        count_bytes=count_bytes,
        count_filename=counts_path.name,
        metadata_bytes=meta_bytes,
        metadata_filename=meta_path.name if meta_path else "",
    )
    results = run_all(ctx)

    expected_ids = set(EXPECTED_FAILS.get(dataset, []))
    n_expected = n_warning_fp = n_error_fp = 0
    any_nonpass = False

    for r in results:
        if r.status == "PASS":
            if verbose:
                print(f"    PASS  {r.rule_id}")
            continue
        if r.status == "SKIP":
            if verbose:
                print(f"    SKIP  {r.rule_id}  {r.message[:60]}")
            continue

        # Non-PASS, non-SKIP result
        any_nonpass = True
        if r.rule_id in expected_ids:
            label = "EXPECTED"
            n_expected += 1
        elif r.severity in ("WARNING",):
            label = "WARNING_FP?" if not strict else "FALSE_POS"
            n_warning_fp += 1
        else:
            label = "FALSE_POS"
            n_error_fp += 1

        print(
            f"    [{label}] {r.rule_id} ({r.severity}) — "
            f"{r.message[:80]}"
        )
        if r.rule_id in RULE_THRESHOLDS:
            print(f"      Threshold: {RULE_THRESHOLDS[r.rule_id]}")
        if r.affected_items:
            for item in r.affected_items[:3]:
                print(f"      affected: {item}")

    if not any_nonpass:
        print(f"    ✅ All rules PASS or SKIP — 0 unexpected results")

    return n_expected, n_warning_fp, n_error_fp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat WARNING-severity results as false positives.",
    )
    parser.add_argument("--dataset", "-d", metavar="NAME")
    args = parser.parse_args()

    if args.dataset:
        datasets = [args.dataset]
    else:
        datasets = sorted(d.name for d in REAL_DIR.iterdir() if d.is_dir())

    print("=" * 65)
    print("BioFlowValidator — Clean Dataset Smoke Test")
    print("=" * 65)
    print(f"Datasets: {datasets}")
    print(f"Strict mode: {args.strict}")
    print()

    total_expected = total_warn = total_error_fp = 0

    for ds in datasets:
        print(f"\n{ds}:")
        n_exp, n_warn, n_err = run_clean_dataset(
            ds, verbose=args.verbose, strict=args.strict
        )
        total_expected += n_exp
        total_warn += n_warn
        total_error_fp += n_err

    print()
    print("=" * 65)
    print("Summary")
    print("=" * 65)
    print(f"  Intentional expected fails : {total_expected}")
    print(f"  WARNING-severity anomalies : {total_warn}")
    print(f"  ERROR false positives       : {total_error_fp}")
    print()

    if total_error_fp == 0 and (not args.strict or total_warn == 0):
        print("✅ PASS — 0 unexpected false positives on clean datasets.")
        print()
        if total_warn > 0 and not args.strict:
            print(
                f"ℹ️  {total_warn} WARNING(s) were detected. Review them above.")
            print(
                "   These may be rule thresholds calibrated too tightly for real "
                "data.\n"
                "   Re-run with --strict to treat them as failures."
            )
        sys.exit(0)
    else:
        fp_count = total_error_fp + (total_warn if args.strict else 0)
        print(f"❌ FAIL — {fp_count} false positive(s) detected on clean datasets.")
        print()
        print(
            "Action required: investigate the rules flagged above.\n"
            "  - If the rule threshold is too tight for real data, raise it.\n"
            "  - If the dataset is genuinely unusual, add it to EXPECTED_FAILS\n"
            "    in this script with a comment explaining why.\n"
            "  - Document any accepted FP in the paper (Table S1) with rationale."
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
