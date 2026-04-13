# GSE37704 — Human HeLa HOXA1 Knockdown

## Provenance

| Field | Value |
|---|---|
| GEO Accession | GSE37704 |
| Paper DOI | https://doi.org/10.1038/nbt.2601 |
| Organism | *Homo sapiens* (GRCh38) |
| Gene ID format | Ensembl (ENSG) |
| Experimental design | HOXA1 knockdown vs control, 3 biological replicates each |
| Replicates | 3 per condition |

## How to Obtain Real Data

1. Navigate to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37704
2. Download supplementary file `GSE37704_FPKM_final.txt.gz` or the raw count matrix.
3. For raw counts, download from the SRA using `sra-tools` or `fasterq-dump` with the SRR accessions listed in the series matrix.
4. Alternatively use the featureCounts output from the recount2 database:
   ```
   https://jhubiostatistics.shinyapps.io/recount/
   ```

## Metadata Extraction
Download `GSE37704_series_matrix.txt.gz` and extract `!Sample_characteristics_ch1` rows.
Standard columns after cleaning: `sample_id`, `condition` (control / HOXA1KD), `batch`.

## Proxy Data Note
`counts.tsv` and `metadata.tsv` are synthetic proxies. Replace with real GEO data before
publication-quality benchmarking.
