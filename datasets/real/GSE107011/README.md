# GSE107011 — Human PBMC Immune Cell Types (Monaco et al. 2019)

## Overview

| Field | Value |
|---|---|
| GEO accession | GSE107011 |
| Organism | *Homo sapiens* (human) |
| Tissue | Peripheral blood |
| Cell types | 29 sorted immune cell populations |
| Design | 29 cell types × 4 donors = 114 samples |
| Gene count | ~20,000 |
| Gene ID format | Gene symbols (HGNC) |
| Library depth | Variable (~15–50 million reads per sample) |
| Count format | Raw integer counts (file: `GSE107011_Processed_data_RCPCmel.txt.gz`) |

**IMPORTANT**: GEO also deposits `GSE107011_Processed_data_TPM.txt.gz`. That
file contains **TPM-normalised values** and must **not** be used as input to
BioFlowValidator — it will (correctly) trigger `NRM-001` (non-integer counts).
Always use the `RCPCmel` file (raw counts per cell type per million-expressed
libraries).

## Publication

> Monaco G, Lee B, Xu W, Mustafah S, Hwang YY, Carré C, Burdin N, Visan L,
> Ceccarelli M, Kemeny L, Andreev M, Han B, Biswas SK, Fairhurst AM, Nardin A,
> Krishnamurthy P, Poidinger M, Newell EW, Larbi A, et al. (2019). RNA-Seq
> Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human
> Cell Types. *Cell Reports*, 26(6), 1627–1640.e7.
> doi:10.1016/j.celrep.2019.01.041

## Download

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107011/suppl/
    GSE107011_Processed_data_RCPCmel.txt.gz
```

To download and standardise automatically:

```bash
python datasets/fetch_real_data.py --dataset GSE107011
```

The script verifies integer counts, parses `<CellType>_<DonorN>` column names
into `condition` (cell type) and `batch` (donor) columns, and writes
`counts.tsv` + `metadata.tsv`.

## Why a Good Benchmark

1. **Scale stress test**: 114 samples exercises `SMP-005` O(n²) correlation
   guard, `NRM-002` library size variation across radically different cell types,
   and `SMP-004` replicate count (4 donors per cell type).

2. **BIO-005 false-positive stress test**: Granulocytes (neutrophils, basophils)
   have physiologically elevated mitochondrial content. This dataset will fire
   `BIO-005` on those samples — document as an expected false positive due to
   tissue biology, not a data error.

3. **BIO-004 dominant-gene test**: Erythroblasts contain haemoglobin genes
   (HBA1, HBA2, HBB) that can exceed 50% of library counts — another expected
   physiological false positive for `BIO-004`.

4. **BIO-007 design test**: Donor (4 donors) as batch vs. cell type (29 types)
   as condition. Because there are 4 × 29 = 116 combinations in a balanced
   design, Cramér's V should be close to 0 — tests the multi-level case.

5. **Rich metadata**: 29 cell-type labels let you test `BIO-001` (multi-condition),
   `BIO-002` (condition label sanity), and `BIO-003` (cardinality) in realistic
   conditions.

## Expected Validator Behaviour (Clean Dataset)

| Rule | Expected | Notes |
|---|---|---|
| GEN-005 | SKIP | Gene symbols (HGNC), not Ensembl IDs |
| BIO-004 | WARNING (FP?) | Haemoglobin dominance in erythroblasts is physiological |
| BIO-005 | WARNING (FP?) | Elevated MT fraction in granulocytes is physiological |
| BIO-006 | PASS | Most housekeeping symbols present |
| BIO-007 | PASS | Balanced 29-type × 4-donor design — no confounding |
| SMP-005 | PASS | 500-sample guard: 114 < 500 so pairwise correlation runs |
| All other rules | PASS | Standard clean dataset |

Add `BIO-004` and `BIO-005` to `EXPECTED_FAILS` in `validate_clean_datasets.py`
with a comment explaining the physiological rationale before using this dataset
in automated CI.

## Metadata Column Standardisation

| Standard column | Source | Values |
|---|---|---|
| `condition` | Column name prefix (before last `_`) | cell type (e.g. `CD4Tcm`) |
| `batch` | Column name suffix (after last `_`) | donor ID (e.g. `Donor1`) |

## Files

| File | Description |
|---|---|
| `counts.tsv` | **Structural proxy v3** (GEO FTP blocked in CI; run fetch script for real data) |
| `metadata.tsv` | Sample metadata — **Structural proxy v3** |
| `sha256sums.txt` | SHA-256 checksums of proxy files |
