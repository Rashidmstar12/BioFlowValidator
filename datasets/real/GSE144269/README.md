# GSE144269 — Single-Condition Dataset (BIO-001 Edge Case)

## Provenance

| Field | Value |
|---|---|
| GEO Accession | GSE144269 |
| Organism | *Homo sapiens* (GRCh38) |
| Gene ID format | Ensembl (ENSG) |
| Experimental design | All samples labeled as a single condition (e.g., tumor) |

## Why This Dataset
Tests the BIO-001 edge case: a real GEO submission where all samples were labeled
with the same condition name, making DE analysis impossible without metadata correction.
This validates that BIO-001 correctly fires on real data, not just synthetic datasets.

## ⚠️ Important: BIO-001 FAIL is Intentional

The validator **will** report:
- `BIO-001` → **ERROR** (single condition detected)

This is the **intended scientific test**, not a false positive. The dataset is
"pathologically clean" — it has exactly one expected ERROR and that ERROR correctly
identifies a real submission problem.

The benchmark infrastructure accounts for this:
- `validate_clean_datasets.py` lists BIO-001 under `EXPECTED_FAILS["GSE144269"]`
- `benchmark_eval.py` excludes this ERROR from the false-positive count

## How to Obtain Real Data

### Pre-download verification (critical)
Before downloading, verify the GEO record still has a single condition:
1. Navigate to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144269
2. Check `!Sample_characteristics_ch1` rows in the series matrix.
3. All samples must share the same condition label for this test to be valid.
   If the record has been updated to include multiple conditions, choose a
   different single-condition GEO accession.

### Step-by-step
1. Download the supplementary count matrix from GEO.
2. Download `GSE144269_series_matrix.txt.gz` for metadata.
3. Confirm all `!Sample_characteristics_ch1` rows list the same condition.
4. Place count matrix at: `datasets/real/GSE144269/raw_counts_download.txt`
5. Run: `python datasets/fetch_real_data.py --dataset GSE144269`

### Alternative single-condition datasets
Any GEO dataset where all samples share one condition label works. Examples:
- Atlas datasets where `condition = "primary tumor"` for all samples
- Time-series where `condition = "time_0h"` for all (wrong metadata annotation)

## Metadata Standardization
```
sample_id   → GSM accession
condition   → single value (e.g., "tumor") for all samples
```
No batch column required (BIO-007 will SKIP — single condition means
Cramér's V cannot be computed anyway).

## Expected Validator Behaviour (clean data)
- `BIO-001` → **ERROR** (intended — single condition)
- All other rules → **PASS** or **SKIP**

## Proxy Data Status
High-fidelity synthetic proxy (v2): 5,010 ENSG genes, 6 samples all labeled "tumor".
Replace with real data before final publication benchmarking.
