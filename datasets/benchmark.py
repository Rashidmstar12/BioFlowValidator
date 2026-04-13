"""Benchmark runner: validates faulty datasets and checks against ground-truth YAML."""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

from app.engine.parser import parse_files
from app.engine.runner import run_all

DATASETS = Path(__file__).parent / "examples"
GROUND_TRUTH = Path(__file__).parent / "ground_truth"

CASES = [
    {
        "name": "mixed_gene_ids",
        "count": DATASETS / "faulty/mixed_gene_ids.tsv",
        "meta": None,
        "truth": GROUND_TRUTH / "mixed_gene_ids.yaml",
    },
    {
        "name": "sample_mismatch",
        "count": DATASETS / "faulty/sample_mismatch_counts.tsv",
        "meta": DATASETS / "faulty/sample_mismatch_meta.tsv",
        "truth": GROUND_TRUTH / "sample_mismatch.yaml",
    },
    {
        "name": "normalized_not_raw",
        "count": DATASETS / "faulty/normalized_not_raw.tsv",
        "meta": None,
        "truth": GROUND_TRUTH / "normalized_not_raw.yaml",
    },
    {
        "name": "too_few_replicates",
        "count": DATASETS / "faulty/too_few_replicates_counts.tsv",
        "meta": DATASETS / "faulty/too_few_replicates_meta.tsv",
        "truth": GROUND_TRUTH / "too_few_replicates.yaml",
    },
    {
        "name": "version_suffix_ids",
        "count": DATASETS / "faulty/version_suffix_ids.tsv",
        "meta": None,
        "truth": GROUND_TRUTH / "version_suffix_ids.yaml",
    },
    # Valid dataset — expect 0 errors
    {
        "name": "valid_dataset",
        "count": DATASETS / "valid/rnaseq_counts_valid.tsv",
        "meta": DATASETS / "valid/metadata_valid.tsv",
        "truth": None,
        "expect_no_errors": True,
    },
]


def run_benchmark() -> None:
    passed = 0
    failed = 0

    for case in CASES:
        name = case["name"]
        count_bytes = Path(case["count"]).read_bytes()
        meta_bytes = Path(case["meta"]).read_bytes() if case.get("meta") else None

        ctx = parse_files(
            count_bytes=count_bytes,
            count_filename=Path(case["count"]).name,
            metadata_bytes=meta_bytes,
            metadata_filename=Path(case["meta"]).name if case.get("meta") else "",
        )
        results = run_all(ctx)
        failed_rules = {r.rule_id: r for r in results if r.status == "FAIL"}

        if case.get("expect_no_errors"):
            error_count = sum(1 for r in results if r.status == "FAIL" and r.severity == "ERROR")
            if error_count == 0:
                print(f"  ✅ {name}: no errors (as expected)")
                passed += 1
            else:
                print(f"  ❌ {name}: expected 0 errors, got {error_count}")
                for r in results:
                    if r.status == "FAIL" and r.severity == "ERROR":
                        print(f"       {r.rule_id}: {r.message[:80]}")
                failed += 1
            continue

        truth = yaml.safe_load(Path(case["truth"]).read_text())
        expected_fails = truth.get("expected_fails", [])

        case_ok = True
        for exp in expected_fails:
            rule_id = exp["rule_id"]
            if rule_id in failed_rules:
                actual_sev = failed_rules[rule_id].severity
                if actual_sev == exp["severity"]:
                    print(f"  ✅ {name}: {rule_id} correctly flagged as {actual_sev}")
                else:
                    print(f"  ⚠️  {name}: {rule_id} flagged but severity {actual_sev} ≠ expected {exp['severity']}")
                    case_ok = False
            else:
                print(f"  ❌ {name}: expected {rule_id} to FAIL but it did not")
                case_ok = False

        if case_ok:
            passed += 1
        else:
            failed += 1

    print(f"\nBenchmark: {passed} passed, {failed} failed out of {passed+failed} cases")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    run_benchmark()
