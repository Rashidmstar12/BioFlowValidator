# GSE89189 — Mouse Macrophage LPS Stimulation

## Provenance

| Field | Value |
|---|---|
| GEO Accession | GSE89189 |
| Organism | *Mus musculus* (GRCm38 / mm10) |
| Gene ID format | Ensembl (ENSMUSG), ~14,000 genes |
| Experimental design | LPS stimulation vs unstimulated macrophages, 3 reps each |
| Replicates | 3 per condition |
| Count source | GEO supplementary count matrix (featureCounts / HTSeq pipeline) |

## Why This Dataset
Tests the mouse organism branch of the validation engine:
- ENSMUSG prefix → GEN-005 detects `organism = "mouse"`
- BIO-005 uses the `mt-` (lowercase) mitochondrial pattern for mouse
- BIO-006 uses the GRCm38 Ensembl housekeeping ID set

## How to Obtain Real Data

### Step-by-step
1. Navigate to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89189
2. Download the supplementary count matrix file (featureCounts or HTSeq output).
3. Verify it contains raw integers (no decimal points):
   ```bash
   head -5 GSE89189_counts.txt | cut -f2-7
   # All values should be integers (no dots)
   ```
4. Verify ENSMUSG gene ID format:
   ```bash
   awk 'NR>1 {print $1}' GSE89189_counts.txt | grep -c "^ENSMUSG"
   # Should equal total gene count (~14,000)
   ```
5. Verify mt- gene presence (mouse MT genes, critical for BIO-005):
   ```bash
   awk '{print $1}' GSE89189_counts.txt | grep "^mt-"
   # Should return ~13 mitochondrial genes (mt-Co1, mt-Nd1, etc.)
   ```
   Note: If the dataset uses gene symbols (not ENSMUSG), mt- genes will be
   present. If it uses ENSMUSG IDs, BIO-005 will SKIP (no mt-prefixed ENSMUSG
   IDs exist — Ensembl IDs for MT genes look like ENSMUSG00000064341).

6. Place the file at: `datasets/real/GSE89189/raw_counts_download.txt`
7. Run: `python datasets/fetch_real_data.py --dataset GSE89189`

### Download via GEOparse (Python)
```python
import GEOparse
gse = GEOparse.get_GEO("GSE89189", destdir="/tmp/geo")
# Find supplementary file URLs
for key, val in gse.metadata.items():
    if "supplementary" in key.lower():
        print(val)
```

## Metadata Standardization
Extract from `GSE89189_series_matrix.txt.gz`:
```
!Sample_title         → sample_id  (GSM2360001 … GSM2360006)
!Sample_characteristics_ch1 (treatment) → condition  (unstimulated / LPS)
```

Expected metadata columns: `sample_id`, `condition`, `batch` (all batch1 —
this is a single-batch experiment; batch column is needed for F4 fault injection).

## Pre-submission Orientation Check
```python
import pandas as pd
df = pd.read_csv('counts.tsv', sep='\t', index_col=0)
assert df.shape[0] > df.shape[1], f"Wrong orientation: {df.shape}"
assert (df == df.astype(int)).all().all(), "Non-integer counts"
print(f"OK: {df.shape[0]} ENSMUSG genes × {df.shape[1]} GSM samples")
```

## Expected Validator Behaviour (clean data)
All rules PASS or SKIP. GEN-005 should report `organism = "mouse"`.

## Proxy Data Status
High-fidelity synthetic proxy (v2): 5,010 ENSMUSG genes, ~5M reads/sample.
Replace with real GEO supplementary file before final publication benchmarking.

