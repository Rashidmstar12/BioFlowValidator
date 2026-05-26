# Changelog

## [1.0.0] - 2026-05-26

### Added
- Rule dependency gating system: rules downstream of SMP-001 and NRM-001 are automatically skipped when a fatal input error is detected
- Minimum gene-count prerequisite guards for SMP-005 (≥500 genes) and BIO-006 (≥1000 genes)
- 4 previously undocumented rules added to README: FMT-008, GEN-005, BIO-007, BIO-008

### Fixed
- BUG-001: Metadata-dependent rules (SMP-002, SMP-004, BIO-001–003, BIO-007) no longer run on mismatched sample IDs after SMP-001 fails
- BUG-002: NRM rules (NRM-002 through NRM-006) no longer cascade after NRM-001 detects pre-normalized data
- BUG-003: README now correctly states 32 validation rules (previously stated 28)
- BUG-004: SMP-005 and BIO-006 no longer produce false positives on small test matrices
- BUG-005: Rat mitochondrial Ensembl ID mapping added to BIO-005 to prevent false-negative checks on rat matrices
- BUG-006: Title-case gene symbol and rat Ensembl ID support added to BIO-006 housekeeping gene checks
- BUG-007: Species-agnostic Ensembl ID version checks in GEN-003 to support mouse, rat, zebrafish, and chicken
- SEC-001: Added chunked upload streaming and auto-cleanup in upload router to prevent OOM Denial-of-Service vulnerabilities

### Validated Against
- 8 public GEO datasets: airway, GSE37704, GSE52778, GSE60450, GSE89189, GSE96870, GSE107011, GSE144269
- 7 synthetic fault categories per dataset (56 total faulty cases): pre-normalized data, sample mismatch, transposed matrix, mixed gene IDs, low replicates, library size imbalance, batch confounding
- False positive rate: 0% on clean data (counts-only mode)
- True positive rate: 100% on all 56 faulty cases
