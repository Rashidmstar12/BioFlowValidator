# GSE52778 — Human Airway Smooth Muscle (Himes et al. 2014)

## Overview

| Field | Value |
|---|---|
| GEO accession | GSE52778 |
| Organism | *Homo sapiens* (human) |
| Tissue | Airway smooth muscle cells |
| Cell types | 4 primary cell lines (N61311, N052611, N080611, N061011) |
| Design | 4 cell lines × 2 treatments (untreated / dexamethasone) |
| Samples | 8 |
| Gene count | ~33,469 |
| Gene ID format | Ensembl gene IDs (ENSG…) |
| Library depth | ~25–35 million reads per sample |
| Count format | Raw integer counts |

## Publication

> Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker RM,
> Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri RA Jr,
> Tantisira KG, Weiss ST, Lu Q. (2014). RNA-Seq transcriptome profiling
> identifies CRISPLD2 as a glucocorticoid responsive gene that modulates
> cytokine function in airway smooth muscle cells.
> *PLOS ONE*, 9(6), e99625. doi:10.1371/journal.pone.0099625

This dataset is also distributed as the Bioconductor `airway` package
(release 1.7.0+, ExperimentHub EH1006).

## Download

The count matrix is available directly from GEO FTP (no R required):

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778/suppl/
    GSE52778_pairedsamplesdge.tsv.gz
```

To download and standardise automatically:

```bash
python datasets/fetch_real_data.py --dataset GSE52778
```

The script strips Ensembl version suffixes, assigns `condition` (untreated/dex)
and `batch` (cell line) columns, and writes `counts.tsv` + `metadata.tsv`.

## Why a Good Benchmark

1. **One of the most cited small RNA-seq benchmarks**: Used in >1,000 papers and
   all major Bioconductor vignettes (DESeq2, edgeR, limma). Any false positive is
   immediately recognised by reviewers.

2. **Paired design**: The 4-cell-line / 2-treatment structure produces a
   natural near-confounding scenario: cell line is both batch and a biological
   covariate. `BIO-007` should fire only when batch = cell line is perfectly
   confounded with treatment, not in the balanced paired design.

3. **Validates proxy-airway results**: This is the true underlying data for the
   `datasets/real/airway/` proxy. If benchmark recall changes materially after
   substituting this real data, the proxy was not representative.

4. **Human Ensembl IDs**: Exercises `GEN-005` (human organism detection),
   `BIO-005` (MT fraction with human MT- prefix), and `BIO-006` (human
   housekeeping Ensembl IDs).

## Expected Validator Behaviour (Clean Dataset)

| Rule | Expected | Notes |
|---|---|---|
| GEN-005 | PASS (human detected) | >80% ENSG… prefix |
| BIO-005 | PASS | MT fraction well below 30% in cell lines |
| BIO-006 | PASS | Human housekeeping genes present |
| BIO-007 | SKIP | No confounding in balanced design; or SKIP if no `batch` col |
| All other rules | PASS | Standard clean dataset |

## Metadata Column Standardisation

| Standard column | Source | Values |
|---|---|---|
| `condition` | Sample name suffix | `untreated`, `dex` |
| `batch` | Sample name prefix | cell line ID (e.g. `N61311`) |

## Files

| File | Description |
|---|---|
| `counts.tsv` | **Structural proxy v3** (GEO FTP blocked in CI; run fetch script for real data) |
| `metadata.tsv` | **Structural proxy v3** (same) |
| `sha256sums.txt` | SHA-256 checksums of proxy files |
