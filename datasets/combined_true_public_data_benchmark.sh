#!/usr/bin/env bash
# combined_true_public_data_benchmark.sh
#
# End-to-end true-public-data benchmark run for GSE52778, GSE60450, GSE107011.
#
# Prerequisites:
#   1. Real data has been ingested for all 3 datasets (see TRUE_PUBLIC_DATA_RUNBOOK.md)
#   2. pip install -r datasets/environment.txt (scipy>=1.11 required)
#   3. Run from the repository root: bash datasets/combined_true_public_data_benchmark.sh
#
# Outputs:
#   datasets/benchmark_results.json
#   datasets/benchmark_results.md
#   datasets/benchmark_results_<timestamp>.json  (timestamped copy)
#   datasets/fault_severity_results.json
#   datasets/fault_severity_results.md

set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-backend}"

# ── Interpreter check ─────────────────────────────────────────────────────
PY3=$(command -v python3 2>/dev/null || true)
if [ -z "$PY3" ]; then
    echo "ERROR: python3 not found in PATH. Install Python 3 and retry."
    exit 1
fi

echo "================================================================"
echo "  BioFlowValidator — True-Public-Data Benchmark"
echo "  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "  Python: $($PY3 --version 2>&1)"
echo "================================================================"
echo ""

# ── scipy preflight check ─────────────────────────────────────────────────
echo "── Preflight: scipy version ─────────────────────────────────────────"
$PY3 - <<'SCIPY_CHECK'
import sys
try:
    import scipy
except ImportError:
    print("  ERROR: scipy not installed. Run: pip install scipy>=1.11")
    sys.exit(1)
try:
    from packaging.version import Version
    ok = Version(scipy.__version__) >= Version("1.11")
except ImportError:
    # packaging not available — fall back to integer tuple comparison
    parts = scipy.__version__.split(".")
    ok = (int(parts[0]), int(parts[1])) >= (1, 11)
if not ok:
    print(f"  ERROR: scipy {scipy.__version__} found; >= 1.11 required for BIO-007.")
    sys.exit(1)
print(f"  ✅  scipy {scipy.__version__}")
SCIPY_CHECK
echo ""

# ── 0. Verify all 3 datasets are marked true_public_data ──────────────────
echo "── Step 0: Verify benchmark_mode ───────────────────────────────────"
all_ok=true
for ds in GSE52778 GSE60450 GSE107011; do
    mode_file="datasets/real/$ds/benchmark_mode.txt"
    if [ ! -f "$mode_file" ]; then
        echo "  ❌  $ds: benchmark_mode.txt missing"
        all_ok=false
    else
        mode=$(cat "$mode_file" | tr -d '[:space:]')
        if [ "$mode" = "true_public_data" ]; then
            echo "  ✅  $ds: $mode"
        else
            echo "  ⚠️   $ds: $mode  (expected true_public_data)"
            echo "        Run: python3 datasets/fetch_real_data.py --dataset $ds"
            echo "        OR:  python3 datasets/ingest_real_data.py --dataset $ds --counts <file>"
            all_ok=false
        fi
    fi
done

if [ "$all_ok" = false ]; then
    echo ""
    echo "ERROR: One or more datasets are not true_public_data."
    echo "See datasets/TRUE_PUBLIC_DATA_RUNBOOK.md for per-dataset instructions."
    exit 1
fi
echo ""

# ── Stale faulty-variants warning ────────────────────────────────────────
if [ -d "datasets/real_faulty" ] && [ -n "$(ls -A datasets/real_faulty 2>/dev/null)" ]; then
    echo "── Warning: stale real_faulty files ────────────────────────────────"
    echo "  datasets/real_faulty/ already contains files from a previous run."
    echo "  Step 2 will overwrite the 3 target-dataset subdirectories, but"
    echo "  other subdirectories will remain and be included in the benchmark."
    echo "  If you want a clean slate: rm -rf datasets/real_faulty"
    echo ""
fi

# ── 1. Validate clean datasets ────────────────────────────────────────────
echo "── Step 1: Clean dataset validation ────────────────────────────────"
python3 datasets/validate_clean_datasets.py 2>&1 | tee /tmp/bfv_clean_validation.log
echo ""

# ── 2. Regenerate faulty variants for the 3 datasets ──────────────────────
echo "── Step 2: Generate faulty variants ────────────────────────────────"
for ds in GSE52778 GSE60450 GSE107011; do
    echo "  Generating faults for $ds …"
    python3 datasets/generate_real_faults.py --dataset "$ds"
done
echo ""

# ── 3. Run full benchmark ─────────────────────────────────────────────────
echo "── Step 3: Run full benchmark ──────────────────────────────────────"
python3 datasets/benchmark_eval.py \
    --json-out datasets/benchmark_results.json \
    --md-out datasets/benchmark_results.md \
    --strict
echo ""

# ── 4. Run fault severity analysis ───────────────────────────────────────
echo "── Step 4: Fault severity analysis ─────────────────────────────────"
python3 datasets/fault_severity_analysis.py
echo ""

# ── 5. Print summary ─────────────────────────────────────────────────────
echo "── Step 5: Summary ─────────────────────────────────────────────────"
python3 - <<'PYEOF'
import json, sys
try:
    d = json.load(open("datasets/benchmark_results.json"))
except FileNotFoundError:
    print("  benchmark_results.json not found")
    sys.exit(1)

print(f"  Gates passed:   {d.get('gates_passed')}")
print(f"  Clean FP count: {d.get('clean_fp_count')}")
print(f"  Synthetic:      {d['synthetic_summary']['passed']}/{d['synthetic_summary']['total']}")
print(f"  Real faulty:    {d['real_summary']['passed']}/{d['real_summary']['total']}")
print()
print("  Dataset modes:")
for k, v in sorted(d.get("dataset_modes", {}).items()):
    marker = "✅" if v == "true_public_data" else "⚙️ "
    print(f"    {marker} {k}: {v}")
print()
print("  Fault-type recall:")
for ft, m in sorted(d.get("fault_type_metrics", {}).items()):
    status = "✅" if m["recall"] >= 0.9 else "❌"
    print(f"    {status} {ft}: recall={m['recall']:.3f}  (TP={m['tp']} FN={m['fn']})")
PYEOF

# ── Save timestamped copy of results ────────────────────────────────────
TIMESTAMP=$(date -u +%Y-%m-%dT%H%M%SZ)
if [ -f "datasets/benchmark_results.json" ]; then
    cp "datasets/benchmark_results.json" "datasets/benchmark_results_${TIMESTAMP}.json"
    echo ""
    echo "  Timestamped copy: datasets/benchmark_results_${TIMESTAMP}.json"
fi

echo ""
echo "================================================================"
echo "  Outputs written:"
echo "    datasets/benchmark_results.json"
echo "    datasets/benchmark_results_${TIMESTAMP}.json  (timestamped copy)"
echo "    datasets/benchmark_results.md"
echo "    datasets/fault_severity_results.json"
echo "    datasets/fault_severity_results.md"
echo ""
echo "  Next steps:"
echo "    1. Fill in datasets/external_benchmark_GSE52778.md"
echo "    2. Fill in datasets/external_benchmark_GSE60450.md"
echo "    3. Fill in datasets/external_benchmark_GSE107011.md"
echo "    4. Update datasets/benchmark_manifest.yaml checksums and source_type"
echo "    5. Review the publication-readiness checklist in TRUE_PUBLIC_DATA_RUNBOOK.md"
echo "================================================================"

