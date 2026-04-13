# BioFlowValidator — Validation Rules Reference

Each rule has a unique ID, category, severity, description, and remediation suggestion.

## Severity Levels

| Level | Meaning |
|---|---|
| `ERROR` | Critical issue that will cause DE analysis to fail or produce incorrect results. |
| `WARNING` | Issue that may affect results; should be reviewed before proceeding. |
| `INFO` | Informational note; no action required. |

---

## Format Rules (`format`)

| Rule ID | Severity | Description |
|---|---|---|
| `FMT-001` | ERROR | File must be UTF-8 encoded. Latin-1 / BOM encoding detected. |
| `FMT-002` | WARNING | Delimiter detection — reports TSV or CSV as detected. |
| `FMT-003` | ERROR | Header row with sample IDs must be present. |
| `FMT-004` | ERROR | Duplicate sample (column) names are not allowed. |
| `FMT-005` | ERROR | All count matrix values must be numeric (no text, NA, or empty cells). |
| `FMT-006` | ERROR | Raw counts must be non-negative (≥ 0). |
| `FMT-007` | WARNING | Sample/gene names should not contain leading/trailing whitespace or special characters. |

---

## Sample Rules (`sample`)

| Rule ID | Severity | Description |
|---|---|---|
| `SMP-001` | ERROR | Sample IDs in count matrix and metadata must match exactly (case-sensitive). |
| `SMP-002` | WARNING | Sample order should be consistent between count matrix and metadata. |
| `SMP-003` | ERROR | Sample IDs must be unique in both files. |
| `SMP-004` | ERROR/WARNING | Each condition must have ≥ 2 replicates (ERROR < 2; WARNING = 2, recommend 3+). |
| `SMP-005` | WARNING | Near-identical sample profiles detected (possible technical replicate confusion). |

---

## Gene ID Rules (`gene`)

| Rule ID | Severity | Description |
|---|---|---|
| `GEN-001` | ERROR | Gene IDs must use a single consistent namespace (Ensembl, Entrez, or gene symbols). |
| `GEN-002` | ERROR | Gene IDs must be unique (no duplicate row names). |
| `GEN-003` | WARNING | Ensembl gene IDs should not include version suffixes (e.g. `.14`). |
| `GEN-004` | WARNING | Gene IDs should be biologically meaningful (not random strings). |

---

## Normalization Rules (`normalization`)

| Rule ID | Severity | Description |
|---|---|---|
| `NRM-001` | ERROR | Count values should be integers. Non-integer values (>1%) suggest pre-normalised data (RPKM/TPM). |
| `NRM-002` | WARNING | Library size ratio (max/min) > 10× indicates potential quality issues. |
| `NRM-003` | ERROR | Samples with total counts < 1,000 indicate failed libraries. |
| `NRM-004` | WARNING | All-zero gene rows should be filtered before DE analysis. |
| `NRM-005` | WARNING | >50% of genes with median count < 1 suggests very sparse or low-depth data. |
| `NRM-006` | WARNING | Suspiciously identical gene count rows suggest data duplication. |

---

## Biology Rules (`biology`)

| Rule ID | Severity | Description |
|---|---|---|
| `BIO-001` | ERROR | At least two conditions must be present for DE analysis. |
| `BIO-002` | WARNING | Purely numeric condition labels may be misinterpreted as continuous covariates. |
| `BIO-003` | WARNING | A metadata column with unique value per sample is likely a sample ID, not a condition. |
| `BIO-004` | WARNING | A gene that accounts for >50% of any sample's library may indicate a mapping artifact. |
| `BIO-005` | WARNING | High mitochondrial gene fraction (>30%) suggests poor sample quality. |
| `BIO-006` | WARNING | Most common housekeeping genes absent (only applies when gene symbols are used). |

---

## Adding New Rules

1. Create a new class in the appropriate `backend/app/rules/<category>/rules.py`.
2. Inherit from `BaseRule`; set `rule_id`, `category`, `severity`, `description`.
3. Implement `run(self, context: ValidationContext) -> RuleResult`.
4. Register the class in `backend/app/engine/registry.py`.
5. Add unit tests in `backend/tests/unit/rules/`.
