# GSE144269 — Single-Condition Dataset (BIO-001 Edge Case)

## Provenance

| Field | Value |
|---|---|
| GEO Accession | GSE144269 |
| Organism | *Homo sapiens* (GRCh38) |
| Gene ID format | Ensembl (ENSG) |
| Experimental design | All samples labeled as a single condition (tumor) |

## Why This Dataset
Tests the BIO-001 edge case: a real GEO submission where all samples were labeled
with the same condition name, making DE analysis impossible without metadata correction.
This validates that BIO-001 correctly fires on real data, not just synthetic datasets.

## How to Obtain Real Data

1. Navigate to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144269
2. Download supplementary count matrix.
3. Verify that all samples share the same condition label in the metadata.

## Expected Validation Result (Clean Dataset)
- BIO-001 → FAIL/ERROR (single condition)
- All other rules → PASS or SKIP

This is a "pathological clean" dataset: it has exactly one expected ERROR and that
ERROR is the intended scientific test, not a false positive.

## Proxy Data Note
`counts.tsv` and `metadata.tsv` are synthetic proxies. Replace with real GEO data
before publication-quality benchmarking.
