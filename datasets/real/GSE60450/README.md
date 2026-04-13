# GSE60450 — Mouse Mammary Gland (limma/voom paper dataset)

## Overview

| Field | Value |
|---|---|
| GEO accession | GSE60450 |
| Organism | *Mus musculus* (mouse) |
| Tissue | Mammary gland |
| Cell types | Basal epithelial, Luminal epithelial |
| Design | 2 cell types × 3 developmental stages (virgin / pregnant / lactating) |
| Samples | 12 |
| Gene count | ~27,179 |
| Gene ID format | NCBI Entrez gene IDs (not Ensembl) |
| Library depth | ~4–8 million reads per sample |
| Count format | Raw integer counts (featureCounts, mm10 NCBI RefSeq) |

## Publication

> Law CW, Chen Y, Shi W, Smyth GK. (2014). voom: Precision weights unlock linear
> model analysis tools for RNA-seq read counts. *Genome Biology*, 15, R29.
> doi:10.1186/gb-2014-15-2-r29

This dataset is distributed with the Bioconductor `edgeR` and `limma` packages
and used in the official workflow vignettes.

## Download

The raw count matrix is available directly from GEO FTP (no R required):

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/suppl/
    GSE60450_Lactation-GenewiseCounts.txt.gz
```

To download and standardise automatically:

```bash
python datasets/fetch_real_data.py --dataset GSE60450
```

The script strips the `Length` annotation column, assigns cell-type condition
labels from sample names, and writes `counts.tsv` + `metadata.tsv`.

## Why a Good Benchmark

1. **Widely known ground truth**: Used in dozens of published workflows; any
   unexpected validator behaviour is immediately identifiable.

2. **Two-factor design**: Cell type (basal/luminal) × developmental stage
   provides a natural test for `BIO-007` (batch confounding), `SMP-004`
   (replicate count), and `BIO-001` (single condition).

3. **Mouse Entrez IDs**: Tests `GEN-001` (gene ID format) and `GEN-005`
   (organism detection). Entrez IDs will cause `GEN-005` to SKIP (no Ensembl
   prefix), which is the correct expected behaviour — document in benchmarks.

4. **No R dependency**: Counts are directly downloadable as a TSV.

## Expected Validator Behaviour (Clean Dataset)

| Rule | Expected | Notes |
|---|---|---|
| GEN-005 | SKIP | Entrez IDs; no Ensembl prefix for organism detection |
| BIO-005 | PASS or SKIP | MT detection requires Ensembl/symbol IDs; may SKIP |
| BIO-006 | SKIP | Housekeeping Ensembl IDs not present; will SKIP |
| All other rules | PASS | Standard clean dataset |

## Metadata Column Standardisation

| Standard column | Source | Values |
|---|---|---|
| `condition` | Sample name parsing | `basal`, `luminal` |
| (no batch) | No batch information in original design | — |

## Files

| File | Description |
|---|---|
| `counts.tsv` | Raw integer counts (genes × samples) — populated by fetch script |
| `metadata.tsv` | Sample metadata with `condition` column — populated by fetch script |
| `sha256sums.txt` | SHA-256 checksums — updated by fetch script |
