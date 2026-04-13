# GSE89189 — Mouse Macrophage LPS Stimulation

## Provenance

| Field | Value |
|---|---|
| GEO Accession | GSE89189 |
| Organism | *Mus musculus* (GRCm38 / mm10) |
| Gene ID format | Ensembl (ENSMUSG) |
| Experimental design | LPS stimulation vs unstimulated macrophages, 3 reps each |
| Replicates | 3 per condition |

## Why This Dataset
Tests the mouse organism branch of the validation engine:
- Ensembl mouse prefix `ENSMUSG` triggers organism detection (GEN-005)
- Mitochondrial gene pattern is `mt-` (lowercase) not `MT-` (BIO-005 organism-aware fix)
- Housekeeping genes use mouse ENSMUSG IDs (BIO-006 Ensembl mapping)

## How to Obtain Real Data

1. Navigate to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89189
2. Download supplementary count matrix file.
3. Verify matrix orientation (genes × samples) before use.

## Proxy Data Note
`counts.tsv` and `metadata.tsv` are synthetic proxies with ENSMUSG IDs. Replace with
real GEO data before publication-quality benchmarking.
