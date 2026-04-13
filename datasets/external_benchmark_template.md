# BioFlowValidator — External Benchmark Report Template

<!-- Copy this file to datasets/external_benchmark_<DATASET>.md for each true-public-data run. -->
<!-- Fill in all sections.  Delete placeholder comments before submitting. -->

## Dataset

| Field | Value |
|---|---|
| **Dataset name** | <!-- e.g. GSE60450 --> |
| **GEO accession** | <!-- e.g. GSE60450 --> |
| **Source type** | <!-- proxy OR true_public_data --> |
| **Paper** | <!-- Full citation with DOI --> |
| **Organism** | <!-- human / mouse --> |
| **Gene ID type** | <!-- Ensembl ENSG / Entrez / HGNC symbol --> |
| **Shape (genes × samples)** | <!-- e.g. 27179 × 12 --> |
| **Count type verified** | <!-- raw integers / TPM / CPM (if not raw, STOP) --> |
| **Matrix orientation** | <!-- genes × samples (rows = genes) → correct --> |
| **Date ingested** | <!-- ISO 8601: 2026-04-15 --> |
| **Ingested by** | <!-- your name or CI --> |
| **SHA-256 counts** | <!-- from datasets/real/<NAME>/sha256sums.txt --> |
| **SHA-256 metadata** | <!-- from datasets/real/<NAME>/sha256sums.txt --> |

---

## Metadata Standardization

| Column in original file | Mapped to | Notes |
|---|---|---|
| <!-- original column name --> | `condition` | <!-- e.g. "treatment" → condition --> |
| <!-- original column name --> | `batch` | <!-- e.g. "plate" → batch --> |
| <!-- other columns kept --> | `(kept as-is)` | <!-- notes --> |

### Metadata caveats

<!-- Describe any issues encountered during metadata standardization. -->
<!-- Examples: -->
<!-- - batch column derived from sample names, not explicit in original file -->
<!-- - GSE107011 real file uses TPM, not raw counts; raw counts require SRA download -->
<!-- - GSE60450 uses Entrez IDs; GEN-005 produces INFO (organism not determinable) -->

---

## Section 1: Clean Dataset False-Positive Summary

Run: `PYTHONPATH=backend python datasets/validate_clean_datasets.py`

| Rule category | Rules run | FAIL count | Unexpected FAILs (FP) | Expected violations |
|---|---|---|---|---|
| format | 8 | | | |
| sample | 5 | | | |
| gene | 5 | | | |
| normalization | 6 | | | |
| biology | 8 | | | |
| **Total** | **32** | | | |

### True false positives (unexpected FAILs on clean data)

<!-- List any unexpected FAIL results here. For each: -->
<!-- - Rule ID and message -->
<!-- - Why it is a false positive (biological explanation) -->
<!-- - Proposed mitigation (tissue-type override, threshold adjustment, etc.) -->

| Rule ID | Message | FP type | Mitigation |
|---|---|---|---|
| <!-- e.g. BIO-005 --> | <!-- message --> | <!-- true FP / expected biology --> | <!-- proposed fix --> |

### Expected biological warnings

<!-- List FAIL results that are expected given the biology of this dataset. -->
<!-- These are NOT counted as FP in the benchmark. -->

- <!-- e.g. BIO-007 WARNING: near-confounded batch (Cramér's V = 0.85) in paired design -->
- <!-- e.g. BIO-005 WARNING: high MT fraction in granulocytes (expected for this cell type) -->

---

## Section 2: Faulty Variant Benchmark Results

Run:
```bash
python datasets/generate_real_faults.py --dataset <DATASET>
PYTHONPATH=backend python datasets/benchmark_eval.py --strict
```

### Per-rule-category metrics

| Category | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|
| format | | | | | | |
| sample | | | | | | |
| gene | | | | | | |
| normalization | | | | | | |
| biology | | | | | | |

### Per-fault-type recall

| Fault type | Rule tested | TP | FN | Recall | Notes |
|---|---|---|---|---|---|
| sample_label_mismatch | SMP-001 | | | | |
| mixed_gene_ids | GEN-001 | | | | |
| transposed_matrix | FMT-008 | | | | |
| batch_confounding | BIO-007 | | | | <!-- requires scipy>=1.11 --> |
| low_replicates | SMP-004 | | | | |
| library_size_imbalance | NRM-002 | | | | |
| normalized_not_raw | NRM-001 | | | | |

### BIO-007 scipy note

<!-- Confirm scipy is installed: python3 -c "import scipy; print(scipy.__version__)" -->
<!-- Without scipy, batch_confounding recall = 0.000 -->
- scipy version: <!-- e.g. 1.17.1 -->

---

## Section 3: Dataset-Specific Caveats

<!-- Describe any caveats specific to this dataset that affect benchmark interpretation. -->

### Raw-count verification
<!-- Confirm whether the dataset contains true raw counts. -->
<!-- For GSE107011: real GEO file is TPM (not raw counts). State what file was used. -->

### Biological false positives expected on real data

<!-- Predict which rules may produce expected biological FAILs that are not data errors. -->
<!-- Examples: -->
<!-- - BIO-005 (single-gene dominance): erythroblasts in GSE107011 have HBA1/HBA2 > 50% -->
<!-- - BIO-006 (MT fraction): granulocytes have elevated MT fraction -->
<!-- - NRM-002 (library size): diverse cell types in GSE107011 have large natural size variation -->

### BIO-007 small-sample note
<!-- If avg expected cell count < 5 in the batch × condition contingency table: -->
<!-- "BIO-007 appends Cochran 1954 small-N caveat. V = X.XXX; avg expected cell = Y.Y" -->

---

## Section 4: Dataset Suitability Classification

| Role | Suitable? | Rationale |
|---|---|---|
| Main clean benchmark | <!-- Yes / No --> | <!-- why --> |
| Stress test / complexity set | <!-- Yes / No --> | <!-- why --> |
| Edge case | <!-- Yes / No --> | <!-- why --> |

---

## Section 5: Baseline Comparison

| Fault type | FastQC recall | ENCODE recall | PCA/boxplot recall | BioFlowValidator recall |
|---|---|---|---|---|
| sample_label_mismatch | 0.000 | estimated | 0.000 | |
| mixed_gene_ids | 0.000 | 0.000 | 0.000 | |
| transposed_matrix | 0.000 | estimated | 1.000 | |
| batch_confounding | 0.000 | estimated | 1.000 | |
| low_replicates | 0.000 | estimated | 0.000 | |
| library_size_imbalance | ~0.5 | estimated | 1.000 | |
| normalized_not_raw | 0.000 | estimated | 0.000 | |

---

## Reproduction Instructions

```bash
# 1. Clone and set up environment
git clone https://github.com/Rashidmstar12/BioFlowValidator
cd BioFlowValidator
pip install -r datasets/environment.txt

# 2. Place downloaded GEO files (do NOT run this in CI — GEO FTP is blocked)
#    Download manually and place at the path below
python datasets/ingest_real_data.py \
    --dataset <DATASET> \
    --counts /path/to/raw_counts.tsv \
    --metadata /path/to/metadata.csv \
    --notes "Downloaded from GEO on <DATE>"

# 3. Verify ingestion
python datasets/ingest_real_data.py --dataset <DATASET> --verify-only

# 4. Validate clean dataset
PYTHONPATH=backend python datasets/validate_clean_datasets.py

# 5. Generate faulty variants
python datasets/generate_real_faults.py --dataset <DATASET>

# 6. Run benchmark
PYTHONPATH=backend python datasets/benchmark_eval.py --strict \
    --json-out datasets/benchmark_results_<DATASET>.json \
    --md-out datasets/benchmark_results_<DATASET>.md

# 7. Run severity analysis
PYTHONPATH=backend python datasets/fault_severity_analysis.py
```

---

*Report template: `datasets/external_benchmark_template.md`*
*Populate with real results and rename to `datasets/external_benchmark_<DATASET>.md`.*
