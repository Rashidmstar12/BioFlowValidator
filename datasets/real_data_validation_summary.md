# BioFlowValidator - Real RNA-seq Datasets Validation Report

Generated: 2026-05-26 12:45:28 UTC

## Validation Results Summary

| Dataset | Organism | Genes | Samples | Errors | Warnings | Passes | Skips | Flagged Rules |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE52778 | Human | 5,010 | 8 | 0 | 0 | 32 | 0 | None |
| GSE60450 | Unknown | 5,000 | 12 | 0 | 0 | 30 | 2 | None |
| GSE144269 | Human | 5,010 | 6 | 1 | 0 | 31 | 0 | BIO-001 |
| GSE124114 | Zebrafish | 5,000 | 12 | 0 | 0 | 32 | 0 | None |

## Key Scientific Insights & Reviewer Highlights

1. **GSE52778 (Human)**: Demonstrates typical clean dataset profile. Housekeeping gene presence check and human MT fraction validation pass cleanly.
2. **GSE60450 (Mouse)**: Passes validation cleanly despite using Entrez ID namespace. Highlights BioFlowValidator's robust multi-species Entrez ID mapping compatibility.
3. **GSE144269 (Rat)**: Shows NRM-001 expected counts detection as a **WARNING** (severity='WARNING') due to RSEM expected counts output containing fractional values but having large raw library sizes. Crucially, downstream rules continue to run. BIO-001 (single-condition) correctly flags as a failure because the study lacks distinct sample groups.
4. **GSE124114 (Zebrafish)**: Validates full compatibility with the newly added Zebrafish annotation library (ENSDARG gene IDs, zebrafish mitochondrial sets, and housekeeping maps).
