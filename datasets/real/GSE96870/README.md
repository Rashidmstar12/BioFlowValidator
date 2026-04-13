# GSE96870 — Multi-Factor Human RNA-seq (Sex × Time)

## Provenance

| Field | Value |
|---|---|
| GEO Accession | GSE96870 |
| Organism | *Homo sapiens* (GRCh38) |
| Gene ID format | Ensembl (ENSG) |
| Experimental design | Multi-factor: sex (male/female) × time (0h/6h/24h), 2 reps = 12 samples |
| Replicates | 2 per sex × time combination |
| Count source | GEO supplementary count matrix |

## Why This Dataset
Tests multi-factor confounding detection rules:
- `batch` column with balanced assignment → BIO-007 correctly PASSES (Cramér's V < 0.7)
- Multiple condition axes → tests BIO-003 (metadata cardinality)
- 12 samples → exercises SMP-004 (replicate count) for each condition group

## ⚠️ Critical for BIO-007 Validation
The clean dataset **must not** have a batch column that is confounded with the
`condition` column. After downloading, compute Cramér's V:

```python
import pandas as pd
from scipy.stats import chi2_contingency
meta = pd.read_csv('metadata.tsv', sep='\t', index_col=0)
ct = pd.crosstab(meta['condition'], meta['batch'])
chi2, _, _, _ = chi2_contingency(ct)
v = (chi2 / (ct.values.sum() * (min(ct.shape) - 1))) ** 0.5
print(f"Cramér's V(condition, batch) = {v:.4f}")
# Must be < 0.7 for BIO-007 to not fire on clean data
```

## How to Obtain Real Data

### Step-by-step
1. Navigate to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96870
2. Download the supplementary count matrix file.
3. Download `GSE96870_series_matrix.txt.gz` for metadata.
4. Extract sex and time from `!Sample_characteristics_ch1` rows.
5. Create batch column from sequencing lane or date (from metadata).
6. Verify Cramér's V < 0.7 between batch and condition.
7. Place count matrix at: `datasets/real/GSE96870/raw_counts_download.txt`
8. Run: `python datasets/fetch_real_data.py --dataset GSE96870`

### Alternative: GSE60450 (simpler design)
If GSE96870 metadata proves difficult to parse, **GSE60450** (lactation mouse study,
limma vignette) is structurally equivalent:
- 2 cell types × 2 conditions, 12 samples, confirmed raw counts
- Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450

## Metadata Standardization
After extraction:
```
sample_id   → GSM or SRR accession
condition   → time (0h / 6h / 24h)
batch       → sequencing lane or run date (must be non-confounded with condition)
sex         → covariate column (optional; can be included as extra column)
```

## Expected Validator Behaviour (clean data)
All rules PASS or SKIP. BIO-007 should PASS (V < 0.7 on balanced design).

## Proxy Data Status
High-fidelity synthetic proxy (v2): 5,010 ENSG genes, 12 samples, 3 conditions,
balanced batch assignment. Replace with real data before final publication benchmarking.

