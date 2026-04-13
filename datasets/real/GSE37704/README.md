# GSE37704 — Human HeLa HOXA1 Knockdown

## Provenance

| Field | Value |
|---|---|
| GEO Accession | GSE37704 |
| SRA Accession | SRP012682 |
| Paper DOI | https://doi.org/10.1038/nbt.2601 |
| Organism | *Homo sapiens* (GRCh38) |
| Gene ID format | Ensembl (ENSG) via recount2/recount3 |
| Experimental design | HOXA1 knockdown vs control, 3 biological replicates each |
| Replicates | 3 per condition |
| Raw count source | recount2 (SRP012682) — GEO deposit is FPKM only |

## ⚠️ Critical Sourcing Note

The original GEO deposit (`GSE37704_FPKM_final.txt.gz`) contains **FPKM-normalized
floating-point values**, not raw counts. If submitted directly to BioFlowValidator,
NRM-001 will correctly FAIL (non-integer counts). This is expected behaviour, not a
false positive — it reflects the real problem with the original submission.

For raw counts, use **recount2** (Collado-Torres et al. 2017, Nature Biotechnology):

## How to Obtain Raw Counts

### Option A — recount2 via R (recommended)
```r
BiocManager::install("recount")
library(recount)

# Download Rail-RNA re-aligned counts for SRP012682 (= GSE37704)
download_study("SRP012682")
load("SRP012682/rse_gene.Rdata")

counts <- assay(rse_gene, "counts")      # Ensembl IDs, integer matrix
meta   <- as.data.frame(colData(rse_gene))

# Standardize metadata columns
samples <- colnames(counts)
# Extract condition from metadata (HOXA1KD vs control)
# The title field contains "siLuc_rep1" (control) or "siHOXA1_rep1" (KD)
condition <- ifelse(grepl("siHOXA1", meta$title), "HOXA1KD", "control")

meta_df <- data.frame(
  sample_id = samples,
  condition = condition,
  batch     = paste0("batch", seq_len(length(samples))),  # interleaved
  stringsAsFactors = FALSE
)

# Write
counts_df <- as.data.frame(counts)
counts_df <- cbind(gene_id = rownames(counts_df), counts_df)
write.table(counts_df, "counts.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(meta_df,   "metadata.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

### Option B — recount3 (newer API)
```r
BiocManager::install("recount3")
library(recount3)
proj <- create_rse_manual(
  project="SRP012682", project_home="data_sources/sra",
  organism="human", annotation="gencode_v26", type="gene"
)
counts <- assay(proj, "raw_counts")
```

### Option C — Manual SRA download
SRR accessions: SRR493366, SRR493367, SRR493368 (control) + SRR493369, SRR493370, SRR493371 (HOXA1KD)
```bash
fasterq-dump SRR493366 SRR493367 SRR493368 SRR493369 SRR493370 SRR493371
# Align with STAR, quantify with featureCounts against GRCh38 Ensembl annotation
```

## Standardization

After download:
```bash
python -c "
import pandas as pd
df = pd.read_csv('counts.tsv', sep='\t', index_col=0)
# Strip Ensembl version suffixes
df.index = df.index.str.replace(r'\.\d+$', '', regex=True)
assert (df == df.astype(int)).all().all(), 'Non-integer: not raw counts'
print(f'OK: {df.shape[0]} genes × {df.shape[1]} samples')
"
sha256sum counts.tsv metadata.tsv > sha256sums.txt
```

## Batch Column
No natural batch column exists (3 biological replicates in a single lane).
The clean dataset uses interleaved batch labels (batch1/batch2/batch3) so
BIO-007 (Cramér's V) does not fire. The F4 fault injector overrides this.

## Expected Validator Behaviour (clean data)
All rules PASS or SKIP on raw-count data. NRM-001 will FAIL on the FPKM file.

## Proxy Data Status
High-fidelity synthetic proxy (v2): 5,010 genes, ~5M reads/sample, non-arithmetic
ENSG IDs. Replace with recount2 raw counts before final publication benchmarking.

