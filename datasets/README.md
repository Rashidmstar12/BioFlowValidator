# BioFlowValidator — Example Datasets

## Valid Datasets (`examples/valid/`)

| File | Description |
|---|---|
| `rnaseq_counts_valid.tsv` | 500 Ensembl gene IDs × 6 samples (3 control + 3 treated), integer raw counts, NB-distributed |
| `metadata_valid.tsv` | Matching metadata with `condition` and `batch` columns |
| `rnaseq_counts_large.tsv` | 500 genes × 24 samples (6 conditions × 4 replicates) |
| `metadata_large.tsv` | Matching metadata for large dataset |

## Faulty Datasets (`examples/faulty/`)

Each file has exactly one injected fault for targeted unit-level testing.

| File(s) | Injected Fault | Rule triggered |
|---|---|---|
| `mixed_gene_ids.tsv` | First 10 rows use HGNC gene symbols; rest use Ensembl IDs | `GEN-001` GeneIDFormatRule |
| `sample_mismatch_counts.tsv` + `sample_mismatch_meta.tsv` | Matrix has 6 samples; metadata has only 4 | `SMP-001` SampleMatchRule |
| `normalized_not_raw.tsv` | TPM-like float values (0.1–500.0) instead of raw integer counts | `NRM-001` NonIntegerCountRule |
| `too_few_replicates_counts.tsv` + `too_few_replicates_meta.tsv` | 1 replicate per condition | `SMP-004` ReplicateCountRule |
| `duplicate_samples.tsv` | Column `ctrl_1` appears twice | `FMT-004` / `SMP-003` |
| `single_condition_counts.tsv` + `single_condition_meta.tsv` | All samples labelled `control` (no second condition) | `BIO-001` SingleConditionRule |
| `version_suffix_ids.tsv` | Ensembl IDs have version suffixes e.g. `ENSG00000000001.7` | `GEN-003` VersionSuffixRule |
| `high_zero_rate.tsv` | ≈65% of gene rows are all-zero | `NRM-004` AllZeroGeneRule |
| `library_size_outlier.tsv` | `ctrl_1` has 50× more counts than other samples | `NRM-002` LibrarySizeRule |
| `corrupt_encoding.csv` | File encoded in Latin-1 with non-ASCII characters | `FMT-001` EncodingRule |

## Perturbation Methodology

All faulty datasets were generated programmatically from the valid dataset (`rnaseq_counts_valid.tsv`):

1. Start from a synthetic NB-distributed count matrix (realistic shape).
2. Apply a single targeted perturbation per file.
3. Document the perturbation type and expected rule trigger.

The generation script is `datasets/generate_datasets.py`.

## Ground-Truth Expected Outcomes

Ground-truth YAML files are located in `datasets/ground_truth/`. Each file lists:
- Which rules should FAIL
- Which rules should PASS
- Which rules should SKIP

These are used by the benchmark runner to compute precision/recall per rule category.
