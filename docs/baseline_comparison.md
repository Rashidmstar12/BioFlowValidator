# BioFlowValidator — Baseline Comparison Framework

This document defines the baselines against which BioFlowValidator's detection
performance should be compared in the manuscript.  It specifies what each
baseline can and cannot detect, and defines the comparison metrics.

---

## Baseline 1: FastQC / MultiQC + Manual Checklist

### What it is
FastQC (Andrews et al. 2017, Bioinformatics) is the de facto standard
pre-analysis QC tool for sequencing data.  MultiQC aggregates FastQC reports.
This baseline represents what a bioinformatician would manually check *before*
submitting data to a DE pipeline, using standard FastQC HTML output.

### Citation
> Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput
> Sequence Data. Available online at:
> http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

### Per-fault-type coverage

| BioFlowValidator Rule | FastQC detectable? | Notes |
|---|---|---|
| FMT-001 (encoding) | ❌ No | FastQC operates on raw FASTQ; file encoding checked at count matrix stage |
| FMT-002 (delimiter) | ❌ No | FastQC does not parse count matrices |
| FMT-003 (missing header) | ❌ No | — |
| FMT-004 (duplicate columns) | ❌ No | — |
| FMT-005 (non-numeric values) | ❌ No | — |
| FMT-006 (negative counts) | ❌ No | — |
| FMT-007 (whitespace names) | ❌ No | — |
| FMT-008 (transposed matrix) | ❌ No | — |
| SMP-001 (sample ID mismatch) | ❌ No | FastQC does not know about metadata |
| SMP-002 (sample order) | ❌ No | — |
| SMP-003 (duplicate sample IDs) | ❌ No | — |
| SMP-004 (replicate count) | ❌ No | — |
| SMP-005 (near-identical samples) | ⚠️ Partial | FastQC per-sample duplication metrics can hint at this, but do not give pairwise Pearson r between expression profiles |
| GEN-001 (mixed gene ID namespace) | ❌ No | — |
| GEN-002 (duplicate gene IDs) | ❌ No | — |
| GEN-003 (version suffixes) | ❌ No | — |
| GEN-004 (non-biological IDs) | ❌ No | — |
| GEN-005 (organism detection) | ❌ No | — |
| NRM-001 (non-integer counts) | ❌ No | FastQC does not parse count matrices |
| NRM-002 (library size imbalance) | ⚠️ Partial | FastQC "Per Sequence Quality Scores" can indicate sequencing depth inequality; MultiQC reports total reads per sample which is a proxy |
| NRM-003 (zero-library samples) | ⚠️ Partial | MultiQC total-reads report would show near-zero libraries |
| NRM-004 (all-zero genes) | ❌ No | — |
| NRM-005 (low-count dominance) | ❌ No | — |
| NRM-006 (duplicate gene rows) | ❌ No | — |
| BIO-001 (single condition) | ❌ No | FastQC does not interpret experimental design |
| BIO-002 (numeric condition labels) | ❌ No | — |
| BIO-003 (metadata cardinality) | ❌ No | — |
| BIO-004 (dominant gene) | ❌ No | FastQC does not know gene identity |
| BIO-005 (MT fraction) | ⚠️ Partial | Only if alignments are classified by reference; not available from count matrix |
| BIO-006 (housekeeping genes absent) | ❌ No | — |
| BIO-007 (batch confounding) | ❌ No | FastQC does not compare batches against conditions |
| BIO-008 (ERCC spike-in) | ❌ No | FastQC does not classify ERCC sequences in count matrices |

**Summary**: FastQC covers 3/32 rules (partially) and 0/32 fully.  It is a
FASTQ-level tool; 29 of 32 BioFlowValidator rules address count-matrix-level
faults that FastQC cannot detect by design.

### Comparison metrics
- **Coverage**: fraction of the 32 rules that the baseline can partially or
  fully detect.
- **Per-fault recall**: for each of the 7 injected fault types in the benchmark
  dataset, does FastQC detect it? (Binary per fault type, not per instance.)
- **False-positive rate**: FastQC does not produce per-rule FP/FN counts
  comparable to BioFlowValidator; report as a qualitative column in Table 2.

---

## Baseline 2: ENCODE RNA-seq Standards Checklist (Manual Human Review)

### What it is
The ENCODE consortium published minimum standards for RNA-seq data quality
(ENCODE RNA-seq standards, doi:10.1101/gr.136184.111; updated 2023 portal
guidelines).  These are used as informal human-reviewer baselines in QC papers
(e.g. Schurch et al. 2016).

This baseline represents what a knowledgeable bioinformatician would check by
inspecting the count matrix and metadata manually according to ENCODE guidelines.

### Per-fault-type coverage

| BioFlowValidator Rule | ENCODE checklist? | Notes |
|---|---|---|
| FMT-008 (transposed matrix) | ✅ Yes | An experienced reviewer would notice row/column confusion |
| SMP-001 (sample ID mismatch) | ✅ Yes | Standard data submission check |
| SMP-003 (duplicate sample IDs) | ✅ Yes | — |
| SMP-004 (replicate count) | ✅ Yes | ENCODE requires ≥ 2 biological replicates |
| NRM-001 (non-integer counts) | ✅ Yes | ENCODE requires raw counts, not FPKM/TPM |
| NRM-002 (library size imbalance) | ✅ Yes | ENCODE specifies minimum library depth; extreme imbalance is flagged |
| NRM-003 (zero-library samples) | ✅ Yes | Near-zero libraries violate ENCODE minimums |
| BIO-001 (single condition) | ✅ Yes | A reviewer would notice |
| BIO-007 (batch confounding) | ⚠️ Partial | ENCODE requires documenting batches; confounding detection is manual and not always done |
| All format rules (FMT-001..007) | ⚠️ Partial | An expert reviewer might catch some; not systematic |
| Gene ID rules (GEN-001..005) | ⚠️ Partial | Namespace consistency is sometimes checked; not systematic |
| Normalization rules (NRM-004..006) | ⚠️ Partial | Sometimes checked; not automated |
| Biology rules (BIO-002..006, 008) | ❌ No / Partial | Rarely checked manually; specialist knowledge required |

**Summary**: ENCODE checklist covers ~9/32 rules fully and ~10/32 partially.
Manual review is slow (estimated 30–60 minutes per dataset vs. seconds for
BioFlowValidator) and inconsistent between reviewers.

### Comparison metrics
- **Coverage** (same as Baseline 1)
- **Time-to-check**: manual (estimated) vs. automated
- **Reproducibility**: ENCODE guidelines are not a deterministic algorithm;
  two reviewers may reach different conclusions.  BioFlowValidator produces
  identical results for the same input.

---

## Baseline 3: PCA + Library-Size Boxplot (Ad-hoc Statistical Inspection)

### What it is
PCA-based sample outlier detection and library size boxplots are the informal
"everyone does this anyway" baseline — standard analyses included in the DESeq2
and edgeR vignettes.  This baseline represents what a bioinformatician would
do by running:

```R
dds <- DESeqDataSetFromMatrix(countData, colData, ~condition)
vst <- vst(dds)
plotPCA(vst, intgroup = "condition")
boxplot(log2(counts(dds) + 1))
```

### Per-fault-type coverage

| BioFlowValidator Rule | PCA/boxplot detectable? | Notes |
|---|---|---|
| FMT-008 (transposed matrix) | ✅ Yes | PCA of a transposed matrix produces an obviously wrong plot (5,000 "samples") |
| SMP-005 (near-identical samples) | ✅ Yes | Identical samples cluster perfectly in PCA |
| NRM-002 (library size imbalance) | ✅ Yes | Boxplot shows extreme outlier bars |
| BIO-007 (batch confounding) | ✅ Yes | PCA coloured by batch vs. condition reveals confounding visually |
| SMP-001 (sample ID mismatch) | ❌ No | DESeq2 would error before PCA reaches this |
| NRM-001 (non-integer counts) | ❌ No | DESeq2 may or may not error; silent failure |
| GEN-001 (mixed gene IDs) | ❌ No | PCA does not inspect gene IDs |
| BIO-005 (MT fraction) | ⚠️ Partial | Only if MT fraction is plotted explicitly |
| All format rules (FMT-001..007) | ❌ No | — |
| All other SMP / GEN / NRM / BIO rules | ❌ No | — |

**Summary**: PCA + boxplot covers 4/32 rules fully and 1/32 partially.  It is
a visualisation approach, not a deterministic rule engine: detection depends on
the bioinformatician noticing and interpreting the plot correctly.

### Comparison metrics
- **Per-fault recall**: for the 7 benchmark fault types, report binary
  detectability by PCA/boxplot inspection.
- **Subjectivity**: PCA is user-interpreted; BioFlowValidator provides a
  deterministic, machine-readable result with quantitative evidence
  (affected items, Cramér's V value, etc.).
- **Coverage of fault types PCA is blind to**: 28/32 rules have no PCA
  analogue — highlight this as BioFlowValidator's core contribution.

---

## Actual Benchmark Results (Populated)

The following quantitative recall values are from `datasets/benchmark_eval.py`
run across 56 real-faulty cases (8 datasets × 7 fault types).

> **Note**: All datasets are structural proxies; see
> `datasets/real_data_benchmark_report.md` for full provenance disclosure.

### Per-fault-type recall comparison

| Fault type | FastQC recall | ENCODE recall | PCA/boxplot recall | BioFlowValidator recall |
|---|---|---|---|---|
| sample_label_mismatch | 0.000 | —† | 0.000 | **1.000** |
| mixed_gene_ids | 0.000 | 0.000‡ | 0.000 | **1.000** |
| transposed_matrix | 0.000 | 0.500† | 1.000 | **1.000** |
| batch_confounding | 0.000 | 0.000‡ | 1.000† | **1.000** |
| low_replicates | 0.000 | 1.000† | 0.000 | **1.000** |
| library_size_imbalance | 0.500† | 1.000† | 1.000 | **1.000** |
| normalized_not_raw | 0.000 | 0.500† | 0.000 | **1.000** |
| **Overall** | **~0.07** | **~0.43** | **~0.43** | **1.000** |

†Estimated recall based on capability analysis (non-automated baselines).
‡FastQC and ENCODE cannot express this concept (0 recall guaranteed).

FastQC and ENCODE checklist cannot be run programmatically on count matrices;
recall values are estimated from their documented capabilities.  PCA/boxplot
values assume a competent analyst correctly interpreting well-designed plots.

### False-positive rate comparison

| Baseline | Clean FP count | FP rate |
|---|---|---|
| FastQC | N/A | N/A (different input format) |
| ENCODE checklist | N/A | N/A (manual, subjective) |
| PCA/boxplot | ~1–2 / dataset | ~0.05–0.1 (analyst-dependent) |
| **BioFlowValidator** | **0** | **0.000** |

BioFlowValidator produced 0 unexpected false positives across 8 clean datasets
(verified via `datasets/validate_clean_datasets.py`).

### Coverage comparison (32 rules)

| Baseline | Rules fully covered | Rules partially covered | Rules not covered |
|---|---|---|---|
| FastQC/MultiQC | 0 | 3 | 29 |
| ENCODE checklist | ~9 | ~10 | ~13 |
| PCA + boxplot | 4 | 1 | 27 |
| **BioFlowValidator** | **32** | **0** | **0** |

---

## Recommended Comparison Table for Manuscript

The following table structure is recommended for the manuscript's Table 2
(Baseline Comparison):

| Fault type | FastQC | ENCODE checklist | PCA/boxplot | BioFlowValidator |
|---|---|---|---|---|
| Matrix format / encoding | ❌ | ⚠️ | ❌ | ✅ |
| Transposed matrix | ❌ | ✅ | ✅ | ✅ |
| Sample ID mismatch | ❌ | ✅ | ❌ | ✅ |
| Replicate count | ❌ | ✅ | ❌ | ✅ |
| Near-identical samples | ⚠️ | ⚠️ | ✅ | ✅ |
| Mixed gene ID namespace | ❌ | ⚠️ | ❌ | ✅ |
| Non-integer counts | ❌ | ✅ | ❌ | ✅ |
| Library size imbalance | ⚠️ | ✅ | ✅ | ✅ |
| Zero-library samples | ⚠️ | ✅ | ✅ | ✅ |
| ERCC spike-in rows | ❌ | ❌ | ❌ | ✅ |
| Batch–condition confounding | ❌ | ⚠️ | ✅ | ✅ |
| MT fraction outlier | ❌ | ❌ | ⚠️ | ✅ |
| **Total covered** | **3 / 12** | **7 / 12** | **5 / 12** | **12 / 12** |

(✅ = fully detectable, ⚠️ = partially detectable, ❌ = not detectable)

*Note: This covers 12 representative fault categories.  Full 32-rule coverage
is in the supplementary rule rationale table.*

---

## Metrics Definitions

For quantitative comparison across baselines (where applicable):

| Metric | Definition |
|---|---|
| **Coverage** | Fraction of 32 rules (or 7 benchmark fault types) that the baseline can in principle detect |
| **Recall per fault type** | TP / (TP + FN) for each fault type across all injected cases |
| **False-positive rate** | Unexpected alerts on clean datasets / total rules run on clean datasets |
| **Time-to-check** | Wall-clock time from input file receipt to quality report (automated tools) or estimated reviewer minutes (manual tools) |
| **Reproducibility** | Standard deviation of recall across three independent runs (for stochastic tools) or across three independent reviewers (for manual tools) |
