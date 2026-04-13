#!/usr/bin/env python3
"""Phase 8 benchmark evaluation script.

Computes TP/FP/FN/precision/recall/F1 per rule and per rule category across
all benchmark datasets (synthetic + real faulty variants).

Usage:
    python datasets/benchmark_eval.py [--json-out PATH] [--md-out PATH] [--strict]

Outputs:
  - Markdown table to stdout (and optionally to --md-out)
  - JSON summary to datasets/benchmark_results.json (or --json-out)

Exit code:
  0  All precision/recall gates pass
  1  One or more gates fail (FP on clean data, or recall < threshold)

Gate definitions:
  - Zero false positives on any clean real dataset (datasets/real/)
  - Recall >= 0.9 per fault category across all real_faulty datasets
    (only enforced when --strict flag is set, to allow incremental development)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).parent
BACKEND = ROOT.parent / "backend"
sys.path.insert(0, str(BACKEND))

from app.engine.parser import parse_files
from app.engine.runner import run_all
from app.models.rule_result import RuleResult

SYNTHETIC_DIR = ROOT / "examples"
SYNTHETIC_TRUTH = ROOT / "ground_truth"
REAL_DIR = ROOT / "real"
REAL_FAULTY_DIR = ROOT / "real_faulty"
DEFAULT_JSON_OUT = ROOT / "benchmark_results.json"

RECALL_GATE = 0.9  # minimum acceptable recall per fault category

# ---------------------------------------------------------------------------
# Rule category mapping
# ---------------------------------------------------------------------------

def _category_of(rule_id: str) -> str:
    prefix = rule_id.split("-")[0]
    return {
        "FMT": "format",
        "SMP": "sample",
        "GEN": "gene",
        "NRM": "normalization",
        "BIO": "biology",
    }.get(prefix, "unknown")


# ---------------------------------------------------------------------------
# Dataset loading helpers
# ---------------------------------------------------------------------------

def _run_validator(count_path: Path, meta_path: Path | None) -> list[RuleResult]:
    count_bytes = count_path.read_bytes()
    meta_bytes = meta_path.read_bytes() if meta_path and meta_path.exists() else None
    ctx = parse_files(
        count_bytes=count_bytes,
        count_filename=count_path.name,
        metadata_bytes=meta_bytes,
        metadata_filename=meta_path.name if meta_path else "",
    )
    return run_all(ctx)


def _load_truth(yaml_path: Path) -> dict:
    return yaml.safe_load(yaml_path.read_text())


# ---------------------------------------------------------------------------
# TP/FP/FN computation
# ---------------------------------------------------------------------------

def _evaluate_case(
    results: list[RuleResult],
    truth: dict,
    strict_fp: bool = False,
) -> dict[str, Any]:
    """Compare validator output against ground truth YAML.

    Parameters
    ----------
    strict_fp:
        If True, any unexpected FAIL is counted as a false positive.
        Set to False for fault-injected variants where cascade failures are
        expected (e.g. a transposed matrix will cause many rules to fire).
        Set to True only for clean datasets.

    Returns a dict with per-rule TP/FP/FN breakdown and a list of violations.
    """
    failed_rules = {r.rule_id: r for r in results if r.status == "FAIL"}
    passed_rules = {r.rule_id: r for r in results if r.status == "PASS"}
    skipped_rules = {r.rule_id: r for r in results if r.status == "SKIP"}

    expected_fails_spec = truth.get("expected_fails", [])
    expected_passes = set(truth.get("expected_passes", []))
    expected_skips = set(truth.get("expected_skips", []))

    expected_fail_ids = {spec["rule_id"] for spec in expected_fails_spec}

    tp_rules, fp_rules, fn_rules = [], [], []
    violations = []

    # Check each expected fail
    for spec in expected_fails_spec:
        rule_id = spec["rule_id"]
        if rule_id in failed_rules:
            actual = failed_rules[rule_id]
            tp_rules.append(rule_id)
            # Severity check
            if actual.severity != spec.get("severity", actual.severity):
                violations.append(
                    f"SEVERITY_MISMATCH: {rule_id} expected {spec['severity']} "
                    f"got {actual.severity}"
                )
            # min_affected_items check
            min_items = spec.get("min_affected_items", 0)
            if min_items and len(actual.affected_items) < min_items:
                violations.append(
                    f"INSUFFICIENT_ITEMS: {rule_id} expected >= {min_items} affected "
                    f"items, got {len(actual.affected_items)}"
                )
        else:
            fn_rules.append(rule_id)
            violations.append(f"FN: {rule_id} expected FAIL but did not fire")

    # Check for unexpected fails (false positives) — only in strict mode
    if strict_fp:
        for rule_id, result in failed_rules.items():
            if rule_id not in expected_fail_ids:
                if rule_id not in expected_skips:
                    fp_rules.append(rule_id)
                    violations.append(
                        f"FP: {rule_id} fired unexpectedly "
                        f"(severity={result.severity}, msg={result.message[:80]!r})"
                    )

    return {
        "tp": tp_rules,
        "fp": fp_rules,
        "fn": fn_rules,
        "violations": violations,
    }


# ---------------------------------------------------------------------------
# Metrics aggregation
# ---------------------------------------------------------------------------

def _compute_metrics(tp: int, fp: int, fn: int) -> dict:
    precision = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
    recall = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else float("nan")
    )
    return {"tp": tp, "fp": fp, "fn": fn,
            "precision": round(precision, 3), "recall": round(recall, 3),
            "f1": round(f1, 3)}


def _get_benchmark_mode(ds_dir: Path) -> str:
    """Read benchmark_mode.txt; return 'proxy' or 'true_public_data'."""
    bm_file = ds_dir / "benchmark_mode.txt"
    if bm_file.exists():
        return bm_file.read_text().strip()
    return "proxy"


# ---------------------------------------------------------------------------
# Clean dataset false positive check
# ---------------------------------------------------------------------------

def check_clean_datasets() -> tuple[int, list[str]]:
    """Run validator on all clean real datasets; return (fp_count, violation_list)."""
    fp_count = 0
    violations = []
    dataset_modes: dict[str, str] = {}

    # Special case: GSE144269 is intentionally single-condition (BIO-001 expected)
    # GSE107011: BIO-004 (haemoglobin dominance in erythroblasts) and BIO-005
    # (elevated MT in granulocytes/neutrophils) are expected physiological signals,
    # not validator false positives. Must stay in sync with EXPECTED_FAILS in
    # validate_clean_datasets.py.
    intentional_fails = {
        "GSE144269": {"BIO-001"},
        "GSE107011": {"BIO-004", "BIO-005"},
    }

    for ds_dir in sorted(REAL_DIR.iterdir()):
        if not ds_dir.is_dir():
            continue
        count_path = ds_dir / "counts.tsv"
        meta_path = ds_dir / "metadata.tsv"
        if not count_path.exists():
            continue

        mode = _get_benchmark_mode(ds_dir)
        dataset_modes[ds_dir.name] = mode

        results = _run_validator(count_path, meta_path if meta_path.exists() else None)
        allowed = intentional_fails.get(ds_dir.name, set())

        for r in results:
            if r.status == "FAIL" and r.severity in ("ERROR", "WARNING"):
                if r.rule_id not in allowed:
                    fp_count += 1
                    violations.append(
                        f"CLEAN_FP [{ds_dir.name}|{mode}]: {r.rule_id} "
                        f"({r.severity}) — {r.message[:100]}"
                    )

    return fp_count, violations, dataset_modes


# ---------------------------------------------------------------------------
# Synthetic benchmark (backward-compatible with original benchmark.py)
# ---------------------------------------------------------------------------

def run_synthetic_benchmark() -> tuple[list[dict], int, int]:
    """Run the original synthetic benchmark cases."""
    cases = []
    faulty_dir = SYNTHETIC_DIR / "faulty"
    valid_dir = SYNTHETIC_DIR / "valid"

    # Discover all YAML files in ground_truth/
    if not SYNTHETIC_TRUTH.exists():
        return cases, 0, 0

    passed = failed = 0
    for truth_file in sorted(SYNTHETIC_TRUTH.glob("*.yaml")):
        name = truth_file.stem
        count_path = faulty_dir / f"{name}.tsv"
        if not count_path.exists():
            # Try with various name patterns
            for candidate in faulty_dir.glob(f"{name}*.tsv"):
                count_path = candidate
                break

        if not count_path.exists():
            continue

        # Look for matching meta file
        meta_path = faulty_dir / f"{name.replace('counts', 'meta')}.tsv"
        meta_path2 = faulty_dir / f"{name}_meta.tsv"
        if not meta_path.exists() and meta_path2.exists():
            meta_path = meta_path2

        try:
            results = _run_validator(
                count_path, meta_path if meta_path.exists() else None
            )
            truth = _load_truth(truth_file)
            eval_result = _evaluate_case(results, truth)
            ok = len(eval_result["fn"]) == 0 and len(eval_result["fp"]) == 0
            if ok:
                passed += 1
            else:
                failed += 1
            cases.append({
                "name": name,
                "source": "synthetic",
                "ok": ok,
                **eval_result,
            })
        except Exception as exc:
            failed += 1
            cases.append({
                "name": name, "source": "synthetic", "ok": False,
                "tp": [], "fp": [], "fn": [],
                "violations": [f"ERROR: {exc}"],
            })
    return cases, passed, failed


# ---------------------------------------------------------------------------
# Real faulty benchmark
# ---------------------------------------------------------------------------

def run_real_benchmark() -> tuple[list[dict], int, int]:
    """Run benchmark against all real faulty variant datasets."""
    cases = []
    passed = failed = 0

    if not REAL_FAULTY_DIR.exists():
        return cases, 0, 0

    for ds_dir in sorted(REAL_FAULTY_DIR.iterdir()):
        if not ds_dir.is_dir():
            continue
        for truth_file in sorted(ds_dir.glob("*.yaml")):
            fault_type = truth_file.stem
            count_path = ds_dir / f"{fault_type}_counts.tsv"
            meta_path = ds_dir / f"{fault_type}_meta.tsv"

            if not count_path.exists():
                continue

            try:
                results = _run_validator(
                    count_path, meta_path if meta_path.exists() else None
                )
                truth = _load_truth(truth_file)
                eval_result = _evaluate_case(results, truth)
                ok = len(eval_result["fn"]) == 0 and len(eval_result["fp"]) == 0
                if ok:
                    passed += 1
                else:
                    failed += 1
                cases.append({
                    "name": f"{ds_dir.name}/{fault_type}",
                    "dataset": ds_dir.name,
                    "fault_type": fault_type,
                    "source": "real_faulty",
                    "ok": ok,
                    **eval_result,
                })
            except Exception as exc:
                failed += 1
                cases.append({
                    "name": f"{ds_dir.name}/{fault_type}",
                    "dataset": ds_dir.name,
                    "fault_type": fault_type,
                    "source": "real_faulty",
                    "ok": False,
                    "tp": [], "fp": [], "fn": [],
                    "violations": [f"ERROR: {exc}"],
                })

    return cases, passed, failed


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _category_metrics(cases: list[dict]) -> dict[str, dict]:
    """Aggregate TP/FP/FN per rule category."""
    cat_tp: dict[str, int] = {}
    cat_fp: dict[str, int] = {}
    cat_fn: dict[str, int] = {}

    for case in cases:
        for rule_id in case.get("tp", []):
            cat = _category_of(rule_id)
            cat_tp[cat] = cat_tp.get(cat, 0) + 1
        for rule_id in case.get("fp", []):
            cat = _category_of(rule_id)
            cat_fp[cat] = cat_fp.get(cat, 0) + 1
        for rule_id in case.get("fn", []):
            cat = _category_of(rule_id)
            cat_fn[cat] = cat_fn.get(cat, 0) + 1

    categories = sorted(set(list(cat_tp) + list(cat_fp) + list(cat_fn)))
    return {
        cat: _compute_metrics(
            cat_tp.get(cat, 0), cat_fp.get(cat, 0), cat_fn.get(cat, 0)
        )
        for cat in categories
    }


def _fault_type_metrics(cases: list[dict]) -> dict[str, dict]:
    """Aggregate recall per fault type (real_faulty only)."""
    ft_tp: dict[str, int] = {}
    ft_fn: dict[str, int] = {}

    for case in cases:
        if case.get("source") != "real_faulty":
            continue
        ft = case.get("fault_type", "unknown")
        ft_tp[ft] = ft_tp.get(ft, 0) + len(case.get("tp", []))
        ft_fn[ft] = ft_fn.get(ft, 0) + len(case.get("fn", []))

    return {
        ft: _compute_metrics(ft_tp.get(ft, 0), 0, ft_fn.get(ft, 0))
        for ft in sorted(set(list(ft_tp) + list(ft_fn)))
    }


def _format_markdown_table(rows: list[list[str]], headers: list[str]) -> str:
    widths = [max(len(str(r[i])) for r in [headers] + rows) for i in range(len(headers))]
    sep = "| " + " | ".join("-" * w for w in widths) + " |"
    header_line = "| " + " | ".join(str(h).ljust(widths[i]) for i, h in enumerate(headers)) + " |"
    lines = [header_line, sep]
    for row in rows:
        lines.append("| " + " | ".join(str(v).ljust(widths[i]) for i, v in enumerate(row)) + " |")
    return "\n".join(lines)


def print_report(
    all_cases: list[dict],
    clean_fp: int,
    clean_violations: list[str],
    dataset_modes: dict[str, str] | None = None,
    strict: bool = False,
) -> tuple[str, dict, bool]:
    """Build Markdown report + JSON summary. Returns (markdown, json_data, passed_gates)."""

    real_cases = [c for c in all_cases if c["source"] == "real_faulty"]
    synth_cases = [c for c in all_cases if c["source"] == "synthetic"]

    cat_metrics = _category_metrics(all_cases)
    ft_metrics = _fault_type_metrics(all_cases)

    lines = ["# BioFlowValidator Benchmark Results", ""]

    # ── Dataset inventory (benchmark mode) ───────────────────────────────
    if dataset_modes:
        lines += ["## Dataset Inventory", ""]
        lines.append("| Dataset | Benchmark Mode |")
        lines.append("|---|---|")
        for ds_name, mode in sorted(dataset_modes.items()):
            marker = "✅ true_public_data" if mode == "true_public_data" else "⚙️  proxy"
            lines.append(f"| {ds_name} | {marker} |")
        has_real = any(m == "true_public_data" for m in dataset_modes.values())
        has_proxy = any(m == "proxy" for m in dataset_modes.values())
        if has_proxy and not has_real:
            lines.append("")
            lines.append(
                "> ⚠️  All datasets are **proxy** (structural synthetic data). "
                "Run `python datasets/ingest_real_data.py` to ingest true GEO data."
            )
        elif has_proxy and has_real:
            lines.append("")
            lines.append(
                "> ⚠️  Mixed benchmark: some datasets are proxies, some are true_public_data. "
                "F1 scores reflect proxy quality for proxy datasets."
            )
        lines.append("")

    # ── Clean dataset FP check ────────────────────────────────────────────
    lines += ["## Clean Dataset False Positive Check", ""]
    if clean_fp == 0:
        lines.append("✅ **0 false positives** on clean real datasets.")
    else:
        lines.append(f"❌ **{clean_fp} false positive(s)** on clean real datasets:")
        for v in clean_violations:
            lines.append(f"  - {v}")
    lines.append("")

    # ── Per-category metrics ──────────────────────────────────────────────
    lines += ["## Per-Rule-Category Metrics (all benchmark cases)", ""]
    cat_rows = [
        [cat, str(m["tp"]), str(m["fp"]), str(m["fn"]),
         f"{m['precision']:.3f}", f"{m['recall']:.3f}", f"{m['f1']:.3f}"]
        for cat, m in sorted(cat_metrics.items())
    ]
    lines.append(_format_markdown_table(
        cat_rows,
        ["Category", "TP", "FP", "FN", "Precision", "Recall", "F1"],
    ))
    lines.append("")

    # ── Per-fault-type recall (real faulty only) ──────────────────────────
    lines += ["## Recall Per Fault Type (real faulty datasets)", ""]
    ft_rows = [
        [ft, str(m["tp"]), str(m["fn"]), f"{m['recall']:.3f}",
         "✅" if (m["recall"] >= RECALL_GATE or m["tp"] + m["fn"] == 0) else "❌"]
        for ft, m in sorted(ft_metrics.items())
    ]
    lines.append(_format_markdown_table(
        ft_rows, ["Fault Type", "TP", "FN", "Recall", "Gate (≥0.9)"]
    ))
    lines.append("")

    # ── Summary counts ────────────────────────────────────────────────────
    synth_pass = sum(1 for c in synth_cases if c["ok"])
    synth_total = len(synth_cases)
    real_pass = sum(1 for c in real_cases if c["ok"])
    real_total = len(real_cases)

    lines += [
        "## Summary",
        f"- Synthetic cases: **{synth_pass}/{synth_total}** passed",
        f"- Real faulty cases: **{real_pass}/{real_total}** passed",
        f"- Clean dataset FP: **{clean_fp}**",
        "",
    ]

    # ── Violations detail ─────────────────────────────────────────────────
    all_violations = []
    for c in all_cases:
        for v in c.get("violations", []):
            all_violations.append(f"[{c['name']}] {v}")

    if all_violations:
        lines += ["## Violations Detail", ""]
        for v in all_violations[:50]:
            lines.append(f"- {v}")
        if len(all_violations) > 50:
            lines.append(f"- … and {len(all_violations) - 50} more")
        lines.append("")

    markdown = "\n".join(lines)

    # ── Gate evaluation ───────────────────────────────────────────────────
    gates_pass = True
    if clean_fp > 0:
        gates_pass = False

    if strict:
        for ft, m in ft_metrics.items():
            if m["tp"] + m["fn"] > 0 and m["recall"] < RECALL_GATE:
                gates_pass = False

    # ── JSON summary ──────────────────────────────────────────────────────
    json_data = {
        "clean_fp_count": clean_fp,
        "clean_violations": clean_violations,
        "dataset_modes": dataset_modes or {},
        "category_metrics": cat_metrics,
        "fault_type_metrics": ft_metrics,
        "synthetic_summary": {"passed": synth_pass, "total": synth_total},
        "real_summary": {"passed": real_pass, "total": real_total},
        "gates_passed": gates_pass,
        "cases": [
            {
                "name": c["name"],
                "source": c["source"],
                "ok": c["ok"],
                "tp": c["tp"],
                "fp": c["fp"],
                "fn": c["fn"],
                "violations": c["violations"],
            }
            for c in all_cases
        ],
    }

    return markdown, json_data, gates_pass


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--json-out", default=str(DEFAULT_JSON_OUT),
        metavar="PATH", help="Path to write JSON results (default: datasets/benchmark_results.json)"
    )
    parser.add_argument(
        "--md-out", default=None,
        metavar="PATH", help="Path to write Markdown report (default: stdout only)"
    )
    parser.add_argument(
        "--strict", action="store_true",
        help="Exit 1 if recall < 0.9 per fault type (in addition to FP gate)"
    )
    parser.add_argument(
        "--skip-clean", action="store_true",
        help="Skip clean dataset FP check (faster iteration)"
    )
    args = parser.parse_args()

    print("Running benchmark evaluation …\n")

    # Clean dataset check
    if args.skip_clean:
        clean_fp, clean_violations, dataset_modes = 0, [], {}
    else:
        print("Checking clean datasets for false positives …")
        clean_fp, clean_violations, dataset_modes = check_clean_datasets()
        print(f"  Clean FP: {clean_fp}\n")

    # Synthetic benchmark
    print("Running synthetic benchmark …")
    synth_cases, synth_pass, synth_fail = run_synthetic_benchmark()
    print(f"  {synth_pass} passed, {synth_fail} failed\n")

    # Real faulty benchmark
    print("Running real faulty benchmark …")
    real_cases, real_pass, real_fail = run_real_benchmark()
    print(f"  {real_pass} passed, {real_fail} failed\n")

    all_cases = synth_cases + real_cases

    markdown, json_data, gates_pass = print_report(
        all_cases, clean_fp, clean_violations, dataset_modes=dataset_modes, strict=args.strict
    )

    print(markdown)

    # Write outputs
    Path(args.json_out).write_text(json.dumps(json_data, indent=2))
    print(f"\nJSON results written to: {args.json_out}")

    if args.md_out:
        Path(args.md_out).write_text(markdown)
        print(f"Markdown report written to: {args.md_out}")

    sys.exit(0 if gates_pass else 1)


if __name__ == "__main__":
    main()
