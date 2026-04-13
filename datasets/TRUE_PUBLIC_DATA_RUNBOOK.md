# BioFlowValidator — True-Public-Data Benchmark Runbook

**Target datasets**: GSE52778, GSE60450, GSE107011  
**Purpose**: Replace the structural proxy datasets with real GEO data and run
the full benchmark outside the CI sandbox (where GEO FTP is blocked).

---

## Prerequisites

```bash
# Clone and enter the repository
git clone https://github.com/Rashidmstar12/BioFlowValidator
cd BioFlowValidator

# Install all dependencies
pip install -r datasets/environment.txt

# Verify scipy (REQUIRED for BIO-007 batch_confounding)
python3 -c "import scipy; print(scipy.__version__)"  # must be >= 1.11

# Verify the backend is importable
PYTHONPATH=backend python3 -c "from app.engine.runner import run_all; print('OK')"
```

---

## Dataset 1 — GSE52778 (Human Airway, Himes et al. 2014)

### What this dataset tests

| Rule | Expected result on clean data |
|---|---|
| FMT-* | PASS (genes × samples, TSV, integer counts) |
| GEN-001 | PASS (Ensembl ENSG IDs) |
| GEN-005 | PASS (human, inferred from ENSG) |
| NRM-001 | PASS (raw integer counts) |
| BIO-007 | ⚠️ POSSIBLE WARNING — paired design: cell line is both batch and biological covariate; Cramér's V depends on exact label assignment |
| BIO-006 | PASS (housekeeping genes present) |

### File(s) to download

| File | URL |
|---|---|
| `GSE52778_pairedsamplesdge.tsv.gz` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778/suppl/GSE52778_pairedsamplesdge.tsv.gz` |

```bash
# Alternative: browse to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778
# → "Supplementary file" section → download GSE52778_pairedsamplesdge.tsv.gz
```

### Where to place the file

```
datasets/raw/GSE52778_pairedsamplesdge.tsv.gz   # (any local path is fine)
```

### Option A — Automated ingestion (recommended)

```bash
# Decompress first
gunzip -k /path/to/GSE52778_pairedsamplesdge.tsv.gz
# or keep .gz and let the script handle it via fetch_real_data.py

# The fetch script downloads + standardises in one step (requires internet):
python datasets/fetch_real_data.py --dataset GSE52778
```

### Option B — Manual ingestion (when fetch_real_data.py cannot reach GEO)

```bash
# 1. Decompress the file
gunzip -k /path/to/GSE52778_pairedsamplesdge.tsv.gz
# → /path/to/GSE52778_pairedsamplesdge.tsv

# 2. Ingest
python datasets/ingest_real_data.py \
    --dataset GSE52778 \
    --counts /path/to/GSE52778_pairedsamplesdge.tsv \
    --notes "Downloaded from GEO FTP $(date +%Y-%m-%d)"
```

### Expected ingestion outputs

```
datasets/real/GSE52778/
    counts.tsv          # ~33,469 genes × 8 samples
    metadata.tsv        # sample_id, condition (untreated/dex), batch (cell line)
    sha256sums.txt      # SHA-256 for counts.tsv and metadata.tsv
    benchmark_mode.txt  # contains: true_public_data
    provenance.json     # source file, timestamp, checksums  (Option B only)
    README.md           # updated with download record (Option A only)
```

### Verify benchmark_mode

```bash
cat datasets/real/GSE52778/benchmark_mode.txt
# expected output: true_public_data
```

### End-to-end run sequence

```bash
# 1. Ingest (Option A or B above)

# 2. Validate clean dataset — expect 0 unexpected FAILs
PYTHONPATH=backend python datasets/validate_clean_datasets.py 2>&1 | grep -E "GSE52778|FAIL|PASS|SKIP"

# 3. Generate faulty variants (7 fault types)
python datasets/generate_real_faults.py --dataset GSE52778

# 4. Run full benchmark
PYTHONPATH=backend python datasets/benchmark_eval.py \
    --json-out datasets/benchmark_results.json \
    --md-out datasets/benchmark_results.md \
    --strict

# 5. Confirm mode in report
grep -A 3 "GSE52778" datasets/benchmark_results.md
```

### Known biological caveats

- **BIO-007**: The paired cell-line × treatment design may produce a WARNING if
  `batch` = cell line and `condition` = treatment labels are assigned such that
  Cramér's V > 0.7. This is a **biological property of the paired design**, not a
  data error. Document the Cramér's V value in the report.
- **BIO-005**: Not expected on clean data (smooth muscle, not granulocytes).
- **NRM-002**: Library sizes are ~25–35M; within-dataset variation is modest.

---

## Dataset 2 — GSE60450 (Mouse Mammary, Law et al. 2014)

### What this dataset tests

| Rule | Expected result on clean data |
|---|---|
| FMT-* | PASS |
| GEN-001 | PASS (Entrez IDs — numeric strings) |
| GEN-005 | INFO/SKIP (Entrez IDs cannot identify organism — expected and correct) |
| NRM-001 | PASS (raw integer counts) |
| BIO-007 | PASS (balanced 2-cell-type × 3-stage design; Cramér's V should be low) |
| SMP-004 | PASS (3 replicates per group) |

### File(s) to download

| File | URL |
|---|---|
| `GSE60450_Lactation-GenewiseCounts.txt.gz` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz` |

```bash
# Alternative: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450
# → "Supplementary file" → GSE60450_Lactation-GenewiseCounts.txt.gz
```

### Where to place the file

```
/path/to/GSE60450_Lactation-GenewiseCounts.txt.gz   # any local path
```

### Option A — Automated ingestion (recommended)

```bash
python datasets/fetch_real_data.py --dataset GSE60450
```

### Option B — Manual ingestion

```bash
# 1. Decompress
gunzip -k /path/to/GSE60450_Lactation-GenewiseCounts.txt.gz
# → GSE60450_Lactation-GenewiseCounts.txt

# 2. Note on "Length" column:
# The raw GEO file has a second column named "Length" (gene annotation lengths).
# ingest_real_data.py detects and drops it automatically and prints a log message.
# You do NOT need to strip it manually before running the ingest command.

# 3. Prepare a metadata CSV with these columns:
#    sample_id, condition, batch
#    sample_id  = column headers from counts file (MCL1.DG, MCL1.DH, ...)
#    condition  = cell type: "basal" or "luminal"
#    batch      = developmental stage: "virgin", "pregnant", "lactating"
#
#    Column → condition/batch mapping (from GEO series matrix):
#    MCL1.DG → basal   / virgin
#    MCL1.DH → basal   / virgin
#    MCL1.DI → basal   / pregnant
#    MCL1.DJ → basal   / pregnant
#    MCL1.DK → basal   / lactating
#    MCL1.DL → basal   / lactating
#    MCL1.LA → luminal / virgin
#    MCL1.LB → luminal / virgin
#    MCL1.LC → luminal / pregnant
#    MCL1.LD → luminal / pregnant
#    MCL1.LE → luminal / lactating
#    MCL1.LF → luminal / lactating

# 4. Ingest (after creating metadata.csv)
python3 datasets/ingest_real_data.py \
    --dataset GSE60450 \
    --counts /path/to/GSE60450_Lactation-GenewiseCounts.txt \
    --metadata /path/to/GSE60450_metadata.csv \
    --notes "Downloaded from GEO FTP $(date +%Y-%m-%d)"
```

### Expected ingestion outputs

```
datasets/real/GSE60450/
    counts.tsv          # 27,179 genes × 12 samples (after removing "Length" column)
    metadata.tsv        # sample_id, condition (basal/luminal), batch (stage)
    sha256sums.txt
    benchmark_mode.txt  # contains: true_public_data
    provenance.json     # (Option B only)
```

### Verify benchmark_mode

```bash
cat datasets/real/GSE60450/benchmark_mode.txt
# expected: true_public_data
```

### End-to-end run sequence

```bash
# 1. Ingest (Option A or B above)

# 2. Validate clean dataset
PYTHONPATH=backend python datasets/validate_clean_datasets.py 2>&1 | grep -E "GSE60450|FAIL|PASS|SKIP"

# 3. Generate faulty variants
python datasets/generate_real_faults.py --dataset GSE60450

# 4. Run benchmark
PYTHONPATH=backend python datasets/benchmark_eval.py \
    --json-out datasets/benchmark_results.json \
    --md-out datasets/benchmark_results.md \
    --strict

# 5. Verify
grep -A 3 "GSE60450" datasets/benchmark_results.md
```

### Known biological caveats

- **GEN-005**: Will produce INFO (not PASS/FAIL) because Entrez IDs cannot be
  mapped to an organism without a database lookup. This is correct behaviour —
  document in the report as "expected INFO, not FP".
- **BIO-007**: Balanced design (each cell type appears in all 3 stages); Cramér's V
  should be near 0. If WARNING fires, check that batch labels reflect developmental
  stage (not a confounded mapping).
- **NRM-002**: Library sizes vary 4–8M per sample; expect PASS.

---

## Dataset 3 — GSE107011 (Human PBMC Immune Cells, Monaco et al. 2019)

### What this dataset tests

| Rule | Expected result on clean data |
|---|---|
| FMT-* | PASS |
| GEN-001 | PASS (HGNC gene symbols — validated as symbol namespace) |
| GEN-005 | PASS/INFO (gene symbols; human organism inferred if ENSG not present) |
| NRM-001 | PASS (RCPCmel file = raw integer counts; TPM file would FAIL) |
| NRM-002 | ⚠️ POSSIBLE WARNING — extreme library size variation across cell types |
| SMP-005 | PASS (500-sample guard; 114 samples is within limit) |
| BIO-007 | ⚠️ WARNING expected — 4 donors × 29 cell types; donor batch may be near-confounded with certain cell type subsets |
| BIO-004 | ⚠️ EXPECTED FAIL on erythroblasts (HBA1/HBA2/HBB > 50% of counts) |
| BIO-005 | ⚠️ EXPECTED FAIL on granulocytes (elevated MT fraction — physiological) |

### ⚠️ Critical file selection warning

GEO provides **two** supplementary files for GSE107011:

| Filename | Content | Use? |
|---|---|---|
| `GSE107011_Processed_data_RCPCmel.txt.gz` | **Raw integer counts** | ✅ YES |
| `GSE107011_Processed_data_TPM.txt.gz` | TPM (normalised floats) | ❌ NO — will correctly trigger NRM-001 |

**Always download `RCPCmel`**, not `TPM`.

### File(s) to download

| File | URL |
|---|---|
| `GSE107011_Processed_data_RCPCmel.txt.gz` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107011/suppl/GSE107011_Processed_data_RCPCmel.txt.gz` |

```bash
# Alternative: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
# → "Supplementary file" → GSE107011_Processed_data_RCPCmel.txt.gz
```

### Option A — Automated ingestion (recommended)

```bash
python datasets/fetch_real_data.py --dataset GSE107011
```

### Option B — Manual ingestion

```bash
# 1. Decompress
gunzip -k /path/to/GSE107011_Processed_data_RCPCmel.txt.gz

# 2. Verify it contains integers (not floats):
head -2 /path/to/GSE107011_Processed_data_RCPCmel.txt | cut -f1-4

# 3. Prepare a metadata CSV (required — metadata is NOT auto-derived from column names).
#    Column names in the RCPCmel file have the format CD4Tcm_Donor1, CD4Tcm_Donor2,
#    etc. but ingest_real_data.py does NOT parse these to produce condition/batch.
#    You must supply a --metadata CSV explicitly if you want BIO-007, SMP-004, and
#    BIO-001 to run (they SKIP without metadata).
#
#    Minimal metadata CSV format:
#      sample_id,condition,batch
#      CD4Tcm_Donor1,CD4Tcm,Donor1
#      CD4Tcm_Donor2,CD4Tcm,Donor2
#      ...
#    Derive the full table from the 29 column headers in the RCPCmel file.

# 4. Ingest
python3 datasets/ingest_real_data.py \
    --dataset GSE107011 \
    --counts /path/to/GSE107011_Processed_data_RCPCmel.txt \
    --metadata /path/to/GSE107011_metadata.csv \
    --notes "Downloaded RCPCmel (raw counts) from GEO FTP $(date +%Y-%m-%d); NOT the TPM file"
```

> **Important — metadata is mandatory for full coverage**: without `--metadata`, rules
> BIO-007 (batch confounding), SMP-004 (condition replication), and BIO-001 (single
> condition) will all SKIP. The `ingest_real_data.py` script does **not** attempt to
> split column names on `_Donor` or any other pattern to synthesise metadata.

### Expected ingestion outputs

```
datasets/real/GSE107011/
    counts.tsv          # ~20,396 genes × 114 samples
    metadata.tsv        # sample_id, condition (cell type), batch (donor)
    sha256sums.txt
    benchmark_mode.txt  # contains: true_public_data
    provenance.json     # (Option B only)
```

### Verify benchmark_mode

```bash
cat datasets/real/GSE107011/benchmark_mode.txt
# expected: true_public_data
```

### End-to-end run sequence

```bash
# 1. Ingest (Option A or B above)

# 2. Validate clean dataset
#    NOTE: BIO-004 and BIO-005 FAILs on granulocyte/erythroblast subsets are EXPECTED
PYTHONPATH=backend python datasets/validate_clean_datasets.py 2>&1 | grep -E "GSE107011|FAIL|PASS"

# 3. Generate faulty variants
python datasets/generate_real_faults.py --dataset GSE107011

# 4. Run benchmark
PYTHONPATH=backend python datasets/benchmark_eval.py \
    --json-out datasets/benchmark_results.json \
    --md-out datasets/benchmark_results.md \
    --strict

# 5. Verify
grep -A 3 "GSE107011" datasets/benchmark_results.md
```

### Known biological caveats

- **NRM-002**: The 114-sample dataset spans radically different cell types
  (e.g. erythroblasts vs. monocytes). Library size variation > 10× is expected
  biologically, so NRM-002 may FAIL on the clean dataset. This is an **expected
  biological false positive** — document the ratio and note it reflects real
  biology, not a sequencing problem.
- **BIO-007 small-N caveat**: With 29 cell types × 4 donors, avg expected cell
  count in the contingency table is < 5 (4 samples per cell type / 4 donors = 1.0).
  BIO-007 will append the Cochran 1954 small-N caveat. Record the appended message
  in the report.
- **BIO-004 / BIO-005**: Expected FAILs in erythroblasts and granulocytes.
  Add these to the `intentional_fails` list in `benchmark_eval.py` if needed, or
  document as biological exceptions in the report.
- **SMP-005 500-sample guard**: 114 samples is within the guard (< 500); the
  pairwise correlation scan will run. Expect PASS since each cell type has only
  4 donors and they are not technical duplicates.

---

## Combined Run — All 3 Datasets

Run this sequence after all 3 datasets have been ingested.

```bash
#!/usr/bin/env bash
# combined_true_public_data_benchmark.sh
# Run from the repository root.

set -e
export PYTHONPATH=backend

echo "=== BioFlowValidator True-Public-Data Benchmark ==="
echo "Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ""

# ── 0. Verify all 3 datasets are true_public_data ──────────────────────────
for ds in GSE52778 GSE60450 GSE107011; do
    mode=$(cat datasets/real/$ds/benchmark_mode.txt 2>/dev/null || echo "MISSING")
    if [ "$mode" != "true_public_data" ]; then
        echo "ERROR: $ds benchmark_mode = '$mode' (expected true_public_data)"
        echo "Run: python datasets/fetch_real_data.py --dataset $ds"
        exit 1
    fi
    echo "✅ $ds: $mode"
done
echo ""

# ── 1. Validate clean datasets ─────────────────────────────────────────────
echo "=== Clean dataset validation ==="
python datasets/validate_clean_datasets.py 2>&1 | tee /tmp/clean_validation.log
echo ""

# ── 2. Regenerate faulty variants ──────────────────────────────────────────
echo "=== Generating faulty variants ==="
for ds in GSE52778 GSE60450 GSE107011; do
    echo "  Generating faults for $ds …"
    python datasets/generate_real_faults.py --dataset $ds
done
echo ""

# ── 3. Run full benchmark ──────────────────────────────────────────────────
echo "=== Running benchmark ==="
python datasets/benchmark_eval.py \
    --json-out datasets/benchmark_results.json \
    --md-out datasets/benchmark_results.md \
    --strict
echo ""

# ── 4. Run severity analysis ───────────────────────────────────────────────
echo "=== Running fault severity analysis ==="
python datasets/fault_severity_analysis.py
echo ""

echo "=== Outputs ==="
echo "  datasets/benchmark_results.json"
echo "  datasets/benchmark_results.md"
echo "  datasets/fault_severity_results.json"
echo "  datasets/fault_severity_results.md"
echo ""
echo "=== Dataset inventory ==="
grep "true_public_data\|proxy" datasets/benchmark_results.md | head -20
```

Save this as `datasets/combined_true_public_data_benchmark.sh` and run:

```bash
chmod +x datasets/combined_true_public_data_benchmark.sh
bash datasets/combined_true_public_data_benchmark.sh 2>&1 | tee /tmp/benchmark_run.log
```

---

## Filling in the Benchmark Report Template

After the combined run, fill in one report per dataset:

```bash
# Copy template
cp datasets/external_benchmark_template.md \
   datasets/external_benchmark_GSE52778.md

# Edit: replace all <!-- placeholders --> with real values from the run
# Key values to fill in:
#   - counts_sha256 / metadata_sha256  → from datasets/real/GSE52778/sha256sums.txt
#   - Source type: true_public_data
#   - Clean-run FP summary → from validate_clean_datasets.py output
#   - Faulty-run TP/FP/FN → from benchmark_results.json, cases where dataset=GSE52778
#   - Category-wise F1 → from benchmark_results.json, category_metrics
```

```bash
# Repeat for GSE60450 and GSE107011
cp datasets/external_benchmark_template.md datasets/external_benchmark_GSE60450.md
cp datasets/external_benchmark_template.md datasets/external_benchmark_GSE107011.md
```

---

## Updating the Benchmark Manifest

After a successful run, update `datasets/benchmark_manifest.yaml` for each dataset:

```yaml
# Example: GSE52778 after true-public-data run
- name: GSE52778
  source_type: true_public_data          # ← change from proxy
  counts_sha256: "<from sha256sums.txt>" # ← fill in actual checksum
  metadata_sha256: "<from sha256sums.txt>"
  last_ingested: "2026-MM-DD"            # ← fill in date
```

```bash
# Get the checksums:
cat datasets/real/GSE52778/sha256sums.txt
cat datasets/real/GSE60450/sha256sums.txt
cat datasets/real/GSE107011/sha256sums.txt
```

---

## Publication-Readiness Checklist

The following must all be true before the project can claim
**true-public-data validation**:

### Required files

- [ ] `datasets/real/GSE52778/benchmark_mode.txt` = `true_public_data`
- [ ] `datasets/real/GSE60450/benchmark_mode.txt` = `true_public_data`
- [ ] `datasets/real/GSE107011/benchmark_mode.txt` = `true_public_data`
- [ ] `datasets/real/GSE52778/provenance.json` or updated `README.md` (with download date)
- [ ] `datasets/real/GSE60450/provenance.json` or updated `README.md`
- [ ] `datasets/real/GSE107011/provenance.json` or updated `README.md`
- [ ] `datasets/real/*/sha256sums.txt` updated with real-data checksums
- [ ] `datasets/benchmark_results.json` — produced from true-public-data run
- [ ] `datasets/benchmark_results.md` — Markdown of the above
- [ ] `datasets/fault_severity_results.json` — severity analysis (re-run after ingestion)
- [ ] `datasets/external_benchmark_GSE52778.md` — filled-in report template
- [ ] `datasets/external_benchmark_GSE60450.md` — filled-in report template
- [ ] `datasets/external_benchmark_GSE107011.md` — filled-in report template

### Required benchmark_results.json fields

- [ ] `dataset_modes` shows `true_public_data` for GSE52778, GSE60450, GSE107011
- [ ] `clean_fp_count` = 0 (or all FAILs are documented biological exceptions)
- [ ] `fault_type_metrics[*].recall` ≥ 0.9 for all fault types
- [ ] `gates_passed` = `true`

### Required manifest updates

- [ ] `benchmark_manifest.yaml`: `source_type: true_public_data` for all 3 datasets
- [ ] `benchmark_manifest.yaml`: `counts_sha256` / `metadata_sha256` filled in
- [ ] `benchmark_manifest.yaml`: `last_run.date` updated
- [ ] `benchmark_manifest.yaml`: `last_run.all_cases_passed: true`

### Report sign-off

- [ ] Each `external_benchmark_*.md` has:
  - [ ] Source type = `true_public_data`
  - [ ] SHA-256 checksums filled in (matches `sha256sums.txt`)
  - [ ] Clean-run FP count (0, or biological exceptions documented)
  - [ ] Faulty-run TP/FP/FN and recall per fault type
  - [ ] Category-wise precision/recall/F1
  - [ ] Dataset-specific caveats section completed
  - [ ] Reproduction instructions verified (another person can reproduce from scratch)

### Scipy verification

- [ ] `python3 -c "import scipy; print(scipy.__version__)"` ≥ 1.11
  (Without this, `batch_confounding` recall = 0.000)

---

## Quick Reference Commands

```bash
# Check all 3 dataset modes at once
for ds in GSE52778 GSE60450 GSE107011; do
    echo "$ds: $(cat datasets/real/$ds/benchmark_mode.txt)"
done

# Get checksums for all 3
for ds in GSE52778 GSE60450 GSE107011; do
    echo "=== $ds ===" && cat datasets/real/$ds/sha256sums.txt
done

# Check benchmark_results.json mode summary
python3 -c "
import json; d=json.load(open('datasets/benchmark_results.json'))
for k,v in d.get('dataset_modes',{}).items():
    print(f'  {k}: {v}')
"

# Check benchmark gates
python3 -c "
import json; d=json.load(open('datasets/benchmark_results.json'))
print('Gates passed:', d.get('gates_passed'))
print('Clean FP:', d.get('clean_fp_count'))
for ft, m in d.get('fault_type_metrics',{}).items():
    status = '✅' if m['recall'] >= 0.9 else '❌'
    print(f'  {status} {ft}: recall={m[\"recall\"]}')
"
```
