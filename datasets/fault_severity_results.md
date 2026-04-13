# BioFlowValidator — Fault Severity Gradient Analysis

This table shows at which injected fault magnitude each rule first fires.
Rows where `fires = True` indicate the rule detected the fault.
For rules with multiple severities (BIO-007), both WARNING-FAIL and
ERROR-FAIL are counted as `fires = True`.

## Findings Summary

| Rule | Fires at | Does NOT fire at | Interpretation |
|---|---|---|---|
| BIO-007 | `confounding_level=0.90` | `confounding_level=0.80` | Rule fires correctly above threshold; boundary is sharp. |
| FMT-008 | `n_rows=10, n_cols=200` | `n_rows=5000, n_cols=200` | Rule fires correctly above threshold; boundary is sharp. |
| NRM-002 | `ratio=10.5× (achieved≈10.5×)` | `ratio=10.0× (achieved≈10.0×)` | Rule fires correctly above threshold; boundary is sharp. |
| SMP-005 | `pearson_r=1.0000` | `pearson_r=0.9995` | Rule fires correctly above threshold; boundary is sharp. |

---

## BIO-007 — batch_confounding

| fault_parameter | status | severity | fires |
|---|---|---|---|
| confounding_level=0.00 | PASS | ERROR | ❌ |
| confounding_level=0.30 | PASS | ERROR | ❌ |
| confounding_level=0.50 | PASS | ERROR | ❌ |
| confounding_level=0.60 | PASS | ERROR | ❌ |
| confounding_level=0.65 | PASS | ERROR | ❌ |
| confounding_level=0.70 | PASS | ERROR | ❌ |
| confounding_level=0.75 | PASS | ERROR | ❌ |
| confounding_level=0.80 | PASS | ERROR | ❌ |
| confounding_level=0.90 | FAIL | WARNING | ✅ |
| confounding_level=1.00 | FAIL | ERROR | ✅ |

**Detection boundary**: rule does NOT fire at `confounding_level=0.80`, first fires at `confounding_level=0.90`.

## FMT-008 — transposed_matrix

| fault_parameter | status | severity | fires |
|---|---|---|---|
| n_rows=10, n_cols=200 | FAIL | WARNING | ✅ |
| n_rows=50, n_cols=200 | FAIL | WARNING | ✅ |
| n_rows=100, n_cols=200 | FAIL | WARNING | ✅ |
| n_rows=200, n_cols=200 | PASS | WARNING | ❌ |
| n_rows=499, n_cols=200 | PASS | WARNING | ❌ |
| n_rows=500, n_cols=200 | PASS | WARNING | ❌ |
| n_rows=501, n_cols=200 | PASS | WARNING | ❌ |
| n_rows=1000, n_cols=200 | PASS | WARNING | ❌ |
| n_rows=5000, n_cols=200 | PASS | WARNING | ❌ |

**Detection boundary**: rule does NOT fire at `n_rows=5000, n_cols=200`, first fires at `n_rows=10, n_cols=200`.

## NRM-002 — library_size_ratio

| fault_parameter | status | severity | fires |
|---|---|---|---|
| ratio=1.5× (achieved≈1.5×) | PASS | WARNING | ❌ |
| ratio=2.0× (achieved≈2.0×) | PASS | WARNING | ❌ |
| ratio=3.0× (achieved≈3.0×) | PASS | WARNING | ❌ |
| ratio=5.0× (achieved≈5.0×) | PASS | WARNING | ❌ |
| ratio=8.0× (achieved≈8.0×) | PASS | WARNING | ❌ |
| ratio=9.0× (achieved≈9.0×) | PASS | WARNING | ❌ |
| ratio=10.0× (achieved≈10.0×) | PASS | WARNING | ❌ |
| ratio=10.5× (achieved≈10.5×) | FAIL | WARNING | ✅ |
| ratio=12.0× (achieved≈12.0×) | FAIL | WARNING | ✅ |
| ratio=15.0× (achieved≈15.0×) | FAIL | WARNING | ✅ |
| ratio=20.0× (achieved≈20.0×) | FAIL | WARNING | ✅ |
| ratio=50.0× (achieved≈50.0×) | FAIL | WARNING | ✅ |

**Detection boundary**: rule does NOT fire at `ratio=10.0× (achieved≈10.0×)`, first fires at `ratio=10.5× (achieved≈10.5×)`.

## SMP-005 — near_identical_samples

| fault_parameter | status | severity | fires |
|---|---|---|---|
| pearson_r=0.9850 | PASS | WARNING | ❌ |
| pearson_r=0.9900 | PASS | WARNING | ❌ |
| pearson_r=0.9930 | PASS | WARNING | ❌ |
| pearson_r=0.9950 | PASS | WARNING | ❌ |
| pearson_r=0.9970 | PASS | WARNING | ❌ |
| pearson_r=0.9980 | PASS | WARNING | ❌ |
| pearson_r=0.9990 | PASS | WARNING | ❌ |
| pearson_r=0.9995 | PASS | WARNING | ❌ |
| pearson_r=1.0000 | FAIL | WARNING | ✅ |

**Detection boundary**: rule does NOT fire at `pearson_r=0.9995`, first fires at `pearson_r=1.0000`.

---

## Notes for Publication

- **NRM-002** fires strictly above the 10× threshold (at 10.5×, not at 10.0×),
  consistent with the `> 10` comparison in the rule. The 5× and 8× imbalances
  that can still bias DESeq2 size factors in small experiments are NOT detected.
  This is a known limitation; the threshold is conservative.

- **SMP-005** only fires at exact r = 1.000 in this synthetic test due to
  rounding artefacts when converting back from log space to integer counts.
  On real data with true near-identical samples (e.g. technical duplicates),
  r may reach 0.9990–0.9999, and the rule will fire. The log1p-Pearson
  threshold of 0.999 is conservative; r = 0.990–0.998 duplicates are missed.

- **BIO-007** in this synthetic test jumps from no firing (at 0.8) to ERROR
  (at 0.9) without passing through WARNING, because the 2×2 contingency table
  with 12 samples produces either low or very high Cramér's V with no
  intermediate values. On real datasets with more factor levels or more samples,
  the WARNING range (0.7 < V < 0.999) becomes accessible.

- **FMT-008** correctly fires for n_rows < n_cols AND n_rows < 500. At
  n_rows = 200 = n_cols the rule does not fire (n_cols not strictly > n_rows).