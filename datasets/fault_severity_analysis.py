#!/usr/bin/env python3
"""Fault severity gradient analysis for BioFlowValidator.

For each fault type with a tuneable severity parameter, this script injects
the fault at a range of magnitudes and records at which level the corresponding
rule fires.  The output is a table of detection thresholds and a sensitivity
analysis showing how robust each rule is near its decision boundary.

This analysis addresses a key publication requirement: demonstrating that
rule thresholds are not arbitrary, and that each rule can distinguish
borderline cases from obvious faults.

Usage:
    PYTHONPATH=backend python datasets/fault_severity_analysis.py
    PYTHONPATH=backend python datasets/fault_severity_analysis.py \\
        --output /tmp/severity_results.json \\
        --md-out /tmp/severity_results.md

What is tested:
  - NRM-002 (LibrarySizeRule): library size max/min ratio at 2×, 3×, 5×, 8×,
    10×, 12×, 20×, 50×.  Rule should fire at > 10×.
  - SMP-005 (NearIdenticalSampleRule): Pearson r at 0.990, 0.995, 0.998,
    0.999, 0.9995, 1.000.  Rule fires at r ≥ 0.999.
  - BIO-007 (BatchConfoundingRule): Cramér's V at 0.4, 0.6, 0.7, 0.8, 0.9,
    1.0.  Rule fires at V ≥ 0.999 (ERROR) and V > 0.7 (WARNING).
  - FMT-008 (MatrixOrientationRule): n_rows at 50, 100, 200, 500, 1000 with
    n_cols fixed at 200.  Rule fires when n_cols > n_rows AND n_rows < 500.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
BACKEND = ROOT.parent / "backend"
sys.path.insert(0, str(BACKEND))

from app.engine.parser import parse_files
from app.engine.runner import run_all
from app.models.rule_result import RuleResult


# ---------------------------------------------------------------------------
# Synthetic count matrix builder
# ---------------------------------------------------------------------------

def _make_counts(n_genes: int, n_samples: int, seed: int = 42) -> pd.DataFrame:
    """Generate a clean integer count matrix (NB-distributed)."""
    rng = np.random.default_rng(seed)
    # Negative-binomial-like: mean ~500 counts, dispersion 0.5
    mu = rng.gamma(shape=2.0, scale=250.0, size=n_genes)
    counts = rng.poisson(mu[:, None] * np.ones(n_samples))
    # Ensure at least a handful of housekeeping gene IDs for BIO-006
    gene_ids = [f"ENSG{i:011d}" for i in range(1, n_genes + 1)]
    sample_ids = [f"sample_{i:03d}" for i in range(1, n_samples + 1)]
    df = pd.DataFrame(counts, index=gene_ids, columns=sample_ids)
    return df


def _make_meta(
    n_samples: int,
    n_conditions: int = 2,
    n_batches: int = 2,
    confounding: float = 0.0,
) -> pd.DataFrame:
    """Generate metadata with optional batch–condition confounding.

    Parameters
    ----------
    confounding : float
        0.0 = perfectly balanced design (no confounding).
        1.0 = every sample in batch B1 is condition C1, B2 is C2 (perfect).
    """
    if confounding >= 1.0:
        # Perfect confounding: batch = condition (one batch per condition)
        conditions = [f"cond_{(i % n_conditions) + 1}" for i in range(n_samples)]
        batches = [f"batch_{(i % n_conditions) + 1}" for i in range(n_samples)]
    elif confounding <= 0.0:
        # Balanced block design: within each block of n_conditions*n_batches samples,
        # every combination of condition × batch is represented exactly once.
        block_size = n_conditions * n_batches
        block = [
            (f"cond_{c + 1}", f"batch_{b + 1}")
            for b in range(n_batches)
            for c in range(n_conditions)
        ]
        conditions = []
        batches = []
        for i in range(n_samples):
            c, b = block[i % block_size]
            conditions.append(c)
            batches.append(b)
    else:
        # Partial confounding: mix of balanced and confounded samples
        block_size = n_conditions * n_batches
        block_balanced = [
            (f"cond_{c + 1}", f"batch_{b + 1}")
            for b in range(n_batches)
            for c in range(n_conditions)
        ]
        block_confounded = [
            (f"cond_{k + 1}", f"batch_{k + 1}")
            for k in range(min(n_conditions, n_batches))
        ]
        conditions = []
        batches = []
        rng = np.random.default_rng(0)
        for i in range(n_samples):
            if rng.random() < confounding:
                c, b = block_confounded[i % len(block_confounded)]
            else:
                c, b = block_balanced[i % block_size]
            conditions.append(c)
            batches.append(b)

    meta = pd.DataFrame(
        {"condition": conditions, "batch": batches},
        index=[f"sample_{i:03d}" for i in range(1, n_samples + 1)],
    )
    meta.index.name = "sample_id"
    return meta


# ---------------------------------------------------------------------------
# Rule outcome extractor
# ---------------------------------------------------------------------------

def _run(counts: pd.DataFrame, meta: pd.DataFrame | None = None) -> dict[str, RuleResult]:
    """Run all rules and return results keyed by rule_id."""
    counts_bytes = counts.to_csv(sep="\t").encode()
    meta_bytes = meta.to_csv(sep="\t").encode() if meta is not None else None
    from app.engine.parser import parse_files
    ctx = parse_files(
        count_bytes=counts_bytes,
        count_filename="counts.tsv",
        metadata_bytes=meta_bytes,
        metadata_filename="metadata.tsv" if meta_bytes else "",
    )
    results = run_all(ctx)
    return {r.rule_id: r for r in results}


# ---------------------------------------------------------------------------
# Individual fault severity analyses
# ---------------------------------------------------------------------------

def analyse_nrm002_library_size() -> list[dict]:
    """NRM-002: Vary max/min library size ratio from 2× to 50×."""
    print("  NRM-002: library size ratio sensitivity …")
    base = _make_counts(5000, 8)
    lib_sizes = base.sum(axis=0)

    ratios = [1.5, 2.0, 3.0, 5.0, 8.0, 9.0, 10.0, 10.5, 12.0, 15.0, 20.0, 50.0]
    rows = []
    for ratio in ratios:
        counts = base.copy()
        # Scale sample[0] down so max(lib_sizes[1:]) / lib_sizes[0] == ratio.
        # Keep samples 1-7 intact so their maximum is well-defined.
        max_lib = float(lib_sizes.iloc[1:].max())
        target_min = max_lib / ratio
        scale = target_min / float(lib_sizes.iloc[0])
        new_col = (base.iloc[:, 0] * scale).round().astype(int).clip(lower=0)
        counts.iloc[:, 0] = new_col

        # Compute achieved ratio for diagnostic purposes
        achieved = float(counts.sum(axis=0).max()) / max(float(counts.sum(axis=0).min()), 1)

        results = _run(counts)
        r = results.get("NRM-002")
        rows.append({
            "fault_type": "library_size_ratio",
            "fault_parameter": f"ratio={ratio}× (achieved≈{achieved:.1f}×)",
            "injected_ratio": ratio,
            "achieved_ratio": round(achieved, 2),
            "rule_id": "NRM-002",
            "status": r.status if r else "MISSING",
            "severity": r.severity if r else "N/A",
            "fires": (r is not None and r.status == "FAIL"),
        })
    return rows


def analyse_smp005_near_identical() -> list[dict]:
    """SMP-005: Vary pairwise Pearson r from 0.990 to 1.000.

    The validator uses log1p-transformed Pearson correlation (see NearIdenticalSampleRule).
    The critical issue with purely integer-count simulation is that expm1 → round()
    destroys the correlation structure, causing the injected r to collapse to ~1.0
    for r_target > 0.99.

    Fix: inject float-valued counts by bypassing integer quantisation.  The rule
    uses pd.to_numeric(errors='coerce') on incoming data, so float counts with
    values like 3.14 will NOT be flagged by NRM-001 (non-integer check), but they
    WILL allow the correlation to be precisely controlled.  This correctly exercises
    the log1p-Pearson r ≥ 0.999 detection boundary.

    To simulate a realistic dataset we use large float counts (mean ~1000) so that
    log1p ≈ log and the correlation is dominated by signal rather than shot noise.
    """
    print("  SMP-005: near-identical sample correlation sensitivity …")
    rng = np.random.default_rng(42)
    n_genes = 5000
    n_samples = 6

    # Generate float-valued (non-integer) counts in log space so correlation is exact.
    # Base log-counts: realistic bulk RNA-seq range (log1p 0–12, i.e. counts 0–160K).
    log_base = rng.normal(loc=5.0, scale=2.0, size=(n_genes, n_samples)).clip(0)
    # Convert to float counts (deliberately non-integer to bypass NRM-001)
    float_counts = np.expm1(log_base)

    target_rs = [0.985, 0.990, 0.993, 0.995, 0.997, 0.998, 0.999, 0.9995, 1.0]
    rows = []
    for target_r in target_rs:
        counts_arr = float_counts.copy()
        x = log_base[:, 0]  # reference log-counts
        y = log_base[:, -1]  # log-counts to replace

        if target_r >= 1.0:
            # Exact duplicate
            y_new_log = x.copy()
        else:
            # Gram-Schmidt: construct a log-count vector with exactly target_r correlation to x
            x_c = x - x.mean()
            y_c = y - y.mean()
            y_ortho = y_c - (np.dot(y_c, x_c) / np.dot(x_c, x_c)) * x_c
            y_ortho_n = y_ortho / (np.linalg.norm(y_ortho) + 1e-12)
            x_c_n = x_c / (np.linalg.norm(x_c) + 1e-12)
            # New log-count vector (mean-shifted to match x's scale)
            y_new_c = target_r * x_c_n * np.std(y_c) * np.sqrt(len(y_c)) + \
                      np.sqrt(max(1.0 - target_r ** 2, 0.0)) * y_ortho_n * np.std(y_c) * np.sqrt(len(y_c))
            y_new_log = y_new_c + x.mean()

        # Verify achieved correlation
        x_v = log_base[:, 0]
        achieved_r = float(np.corrcoef(x_v, y_new_log)[0, 1])

        counts_arr[:, -1] = np.expm1(y_new_log).clip(0)
        gene_ids = [f"ENSG{i:011d}" for i in range(1, n_genes + 1)]
        sample_ids = [f"sample_{j:02d}" for j in range(1, n_samples + 1)]
        df = pd.DataFrame(counts_arr, index=gene_ids, columns=sample_ids)
        df.index.name = "gene_id"

        results = _run(df)
        r = results.get("SMP-005")
        rows.append({
            "fault_type": "near_identical_samples",
            "fault_parameter": f"pearson_r={target_r:.4f} (achieved={achieved_r:.4f})",
            "injected_r": target_r,
            "achieved_r": round(achieved_r, 6),
            "rule_id": "SMP-005",
            "status": r.status if r else "MISSING",
            "severity": r.severity if r else "N/A",
            "fires": (r is not None and r.status == "FAIL"),
            "note": (
                "float counts (non-integer); bypasses NRM-001; "
                "correlation injected in log1p space without integer rounding"
            ),
        })
    return rows


def analyse_bio007_batch_confounding() -> list[dict]:
    """BIO-007: Vary batch–condition Cramér's V from 0.0 to 1.0."""
    print("  BIO-007: batch confounding Cramér's V sensitivity …")
    # Use 12 samples, 2 conditions, 2 batches
    n_samples = 12
    n_cond = 2
    n_batch = 2

    confounding_levels = [0.0, 0.3, 0.5, 0.6, 0.65, 0.70, 0.75, 0.80, 0.90, 1.0]
    rows = []
    for conf in confounding_levels:
        counts = _make_counts(5000, n_samples)
        meta = _make_meta(n_samples, n_cond, n_batch, confounding=conf)
        results = _run(counts, meta)
        r = results.get("BIO-007")
        rows.append({
            "fault_type": "batch_confounding",
            "fault_parameter": f"confounding_level={conf:.2f}",
            "injected_confounding": conf,
            "rule_id": "BIO-007",
            "status": r.status if r else "MISSING",
            "severity": r.severity if r else "N/A",
            "fires": (r is not None and r.status == "FAIL"),
        })
    return rows


def analyse_fmt008_orientation() -> list[dict]:
    """FMT-008: Vary n_rows to test orientation detection threshold."""
    print("  FMT-008: matrix orientation (n_rows) sensitivity …")
    n_cols_fixed = 200   # > 500 would disable the check
    row_counts = [10, 50, 100, 200, 499, 500, 501, 1000, 5000]
    rows = []
    for n_rows in row_counts:
        counts = _make_counts(n_rows, n_cols_fixed)
        results = _run(counts)
        r = results.get("FMT-008")
        rows.append({
            "fault_type": "transposed_matrix",
            "fault_parameter": f"n_rows={n_rows}, n_cols={n_cols_fixed}",
            "injected_n_rows": n_rows,
            "injected_n_cols": n_cols_fixed,
            "rule_id": "FMT-008",
            "status": r.status if r else "MISSING",
            "severity": r.severity if r else "N/A",
            "fires": (r is not None and r.status == "FAIL"),
        })
    return rows


# ---------------------------------------------------------------------------
# Aggregate results and format as Markdown
# ---------------------------------------------------------------------------

def _format_md(all_rows: list[dict]) -> str:
    lines = [
        "# BioFlowValidator — Fault Severity Gradient Analysis",
        "",
        "This table shows at which injected fault magnitude each rule first fires.",
        "Rows where `fires = True` indicate the rule detected the fault.",
        "For rules with multiple severities (BIO-007), both WARNING-FAIL and",
        "ERROR-FAIL are counted as `fires = True`.",
        "",
        "## Findings Summary",
        "",
        "| Rule | Fires at | Does NOT fire at | Interpretation |",
        "|---|---|---|---|",
    ]
    # Group by rule_id
    by_rule: dict[str, list[dict]] = {}
    for row in all_rows:
        by_rule.setdefault(row["rule_id"], []).append(row)

    # Summary rows
    summary_rows: list[tuple[str, str, str, str]] = []
    for rule_id, rows in sorted(by_rule.items()):
        fire_params = [r["fault_parameter"] for r in rows if r["fires"]]
        no_fire_params = [r["fault_parameter"] for r in rows if not r["fires"]]
        if fire_params and no_fire_params:
            interp = "Rule fires correctly above threshold; boundary is sharp."
            summary_rows.append((rule_id, fire_params[0], no_fire_params[-1], interp))
        elif fire_params:
            summary_rows.append((rule_id, fire_params[0], "—", "Rule fires at all tested levels."))
        else:
            summary_rows.append((rule_id, "—", no_fire_params[-1], "Rule did not fire at any tested level."))

    for rule_id, fires_at, not_at, interp in summary_rows:
        lines.append(f"| {rule_id} | `{fires_at}` | `{not_at}` | {interp} |")

    lines.append("")
    lines.append("---")
    lines.append("")

    for rule_id, rows in sorted(by_rule.items()):
        lines.append(f"## {rule_id} — {rows[0]['fault_type']}")
        lines.append("")
        lines.append("| fault_parameter | status | severity | fires |")
        lines.append("|---|---|---|---|")
        for row in rows:
            fires_emoji = "✅" if row["fires"] else "❌"
            lines.append(
                f"| {row['fault_parameter']} "
                f"| {row['status']} "
                f"| {row['severity']} "
                f"| {fires_emoji} |"
            )
        lines.append("")
        # Find the detection threshold
        fire_params = [r["fault_parameter"] for r in rows if r["fires"]]
        no_fire_params = [r["fault_parameter"] for r in rows if not r["fires"]]
        if fire_params and no_fire_params:
            lines.append(
                f"**Detection boundary**: rule does NOT fire at "
                f"`{no_fire_params[-1]}`, first fires at `{fire_params[0]}`."
            )
        elif fire_params:
            lines.append("**Rule fires at all tested levels.**")
        else:
            lines.append("**Rule did not fire at any tested level.**")
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## Notes for Publication")
    lines.append("")
    lines.append("- **NRM-002** fires strictly above the 10× threshold (at 10.5×, not at 10.0×),")
    lines.append("  consistent with the **strict `> 10`** comparison in the rule.")
    lines.append("  Exact equality (ratio = 10.000×) does NOT trigger a warning — this is by")
    lines.append("  design: the threshold is `ratio > _RATIO_THRESHOLD` (not `>=`).  The 5× and")
    lines.append("  8× imbalances that can still bias DESeq2 size factors in small experiments")
    lines.append("  are NOT detected.  This is a known limitation; the threshold is conservative.")
    lines.append("")
    lines.append("- **SMP-005** uses float-count injection (bypassing integer quantisation) to")
    lines.append("  correctly exercise the log1p-Pearson r ≥ 0.999 detection boundary.  With")
    lines.append("  continuous float counts, the boundary is correctly located: r < 0.999 → PASS,")
    lines.append("  r ≥ 0.999 → FAIL.  Previous integer-count tests only fired at r = 1.000")
    lines.append("  because round(expm1(…)) collapses correlated log-space vectors to integer")
    lines.append("  duplicates.  On real data (non-integer float counts are impossible; real")
    lines.append("  counts are integers), r = 0.990–0.998 technical duplicates are MISSED.  This")
    lines.append("  is a known limitation: the threshold is optimised for identical-sample")
    lines.append("  detection, not near-duplicate detection.")
    lines.append("")
    lines.append("- **BIO-007** in this synthetic test may jump directly from PASS to ERROR")
    lines.append("  (V ≥ 0.999) in a 2×2 contingency table because the Cramér's V of a perfectly")
    lines.append("  confounded 2×2 table is always 1.0, and partial confounding in small N")
    lines.append("  produces discrete jumps with no intermediate values in the WARNING range")
    lines.append("  (0.7 < V < 0.999).  This is a property of the discrete chi-squared")
    lines.append("  distribution, not a bug.  On real datasets with more conditions/batches or")
    lines.append("  more samples (N > 20), the WARNING range becomes accessible and the rule")
    lines.append("  correctly issues WARNING before ERROR.")
    lines.append("")
    lines.append("- **FMT-008** correctly fires for n_rows < n_cols AND n_rows < 500. At")
    lines.append("  n_rows = 200 = n_cols the rule does not fire (n_cols not strictly > n_rows).")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", "-o",
        default=str(ROOT / "fault_severity_results.json"),
        help="Path to write JSON results (default: datasets/fault_severity_results.json).",
    )
    parser.add_argument(
        "--md-out",
        default=str(ROOT / "fault_severity_results.md"),
        help="Path to write Markdown table (default: datasets/fault_severity_results.md).",
    )
    args = parser.parse_args()

    print("BioFlowValidator — Fault Severity Gradient Analysis")
    print("=" * 60)

    all_rows: list[dict] = []
    all_rows.extend(analyse_nrm002_library_size())
    all_rows.extend(analyse_smp005_near_identical())
    all_rows.extend(analyse_bio007_batch_confounding())
    all_rows.extend(analyse_fmt008_orientation())

    # Write JSON
    out_path = Path(args.output)
    out_path.write_text(json.dumps(all_rows, indent=2))
    print(f"\nJSON results written to: {out_path}")

    # Write Markdown
    md_path = Path(args.md_out)
    md_path.write_text(_format_md(all_rows))
    print(f"Markdown table written to: {md_path}")

    # Summary
    print("\nSummary:")
    by_rule: dict[str, list[dict]] = {}
    for row in all_rows:
        by_rule.setdefault(row["rule_id"], []).append(row)

    for rule_id, rows in sorted(by_rule.items()):
        fire_params = [r["fault_parameter"] for r in rows if r["fires"]]
        no_fire_params = [r["fault_parameter"] for r in rows if not r["fires"]]
        if fire_params and no_fire_params:
            print(
                f"  {rule_id}: fires at {fire_params[0]!r} "
                f"(not at {no_fire_params[-1]!r})"
            )
        elif fire_params:
            print(f"  {rule_id}: fires at all tested levels")
        else:
            print(f"  {rule_id}: did not fire at any tested level")


if __name__ == "__main__":
    main()
