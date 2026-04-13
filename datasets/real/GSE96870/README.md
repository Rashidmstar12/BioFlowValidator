# GSE96870 — Multi-Factor Human RNA-seq (Sex × Time)

## Provenance

| Field | Value |
|---|---|
| GEO Accession | GSE96870 |
| Organism | *Homo sapiens* (GRCh38) |
| Gene ID format | Ensembl (ENSG) |
| Experimental design | Multi-factor: sex (male/female) × time (0h/6h/24h), balanced |
| Replicates | 2 per sex × time combination = 12 samples total |

## Why This Dataset
Tests multi-factor confounding detection rules:
- `batch` column tests BIO-007 (Cramér's V batch-condition association)
- Multiple condition axes test BIO-003 (metadata cardinality)

## How to Obtain Real Data

1. Navigate to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96870
2. Download supplementary count matrix.
3. Extract metadata: sex and time columns from series matrix file.

## Proxy Data Note
`counts.tsv` and `metadata.tsv` are synthetic proxies. Replace with real GEO data
before publication-quality benchmarking.
