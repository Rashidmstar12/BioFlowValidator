"""Gene ID validation rules."""
from __future__ import annotations

import re

import pandas as pd

from app.models.context import ValidationContext
from app.models.rule_result import RuleResult
from app.rules.base import BaseRule

_ENSEMBL_RE = re.compile(r'^ENSG\d{11}(\.\d+)?$', re.IGNORECASE)
_ENSEMBL_BASE_RE = re.compile(r'^ENSG\d{11}$', re.IGNORECASE)
_ENSEMBL_VERSION_RE = re.compile(r'^ENSG\d{11}\.\d+$', re.IGNORECASE)
_ENTREZ_RE = re.compile(r'^\d+$')
_SYMBOL_RE = re.compile(r'^[A-Z][A-Z0-9\-]{1,}$')


def _classify_id(gene_id: str) -> str:
    g = gene_id.strip()
    if _ENSEMBL_BASE_RE.match(g):
        return "ensembl"
    if _ENSEMBL_VERSION_RE.match(g):
        return "ensembl_versioned"
    if _ENTREZ_RE.match(g):
        return "entrez"
    if _SYMBOL_RE.match(g):
        return "symbol"
    return "unknown"


class GeneIDFormatRule(BaseRule):
    rule_id = "GEN-001"
    category = "gene"
    severity = "ERROR"
    description = "Gene IDs must use a single consistent namespace (Ensembl, Entrez, or gene symbols)."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")

        gene_ids = [str(i) for i in context.count_matrix.index[:2000]]
        if not gene_ids:
            return self._skip("Count matrix has no rows.")

        type_counts: dict[str, int] = {}
        for gid in gene_ids:
            t = _classify_id(gid)
            type_counts[t] = type_counts.get(t, 0) + 1

        # Remove versioned from ensembl count (handled by VersionSuffixRule)
        meaningful = {k: v for k, v in type_counts.items() if k != "unknown"}
        # Collapse ensembl and ensembl_versioned as same namespace for mixing detection
        ensembl_total = meaningful.get("ensembl", 0) + meaningful.get("ensembl_versioned", 0)
        non_ensembl = {k: v for k, v in meaningful.items() if k not in ("ensembl", "ensembl_versioned")}

        namespaces_present = []
        if ensembl_total > 0:
            namespaces_present.append("ensembl")
        namespaces_present.extend(non_ensembl.keys())

        if len(namespaces_present) > 1:
            return self._fail(
                f"Mixed gene ID namespaces detected: {namespaces_present}.",
                affected_items=[f"{k}: {v} IDs" for k, v in type_counts.items()],
                suggestion=(
                    "Use a single gene ID namespace throughout the count matrix. "
                    "Recommended: Ensembl IDs (ENSG...)."
                ),
                details={"namespace_counts": type_counts},
            )

        dominant = max(type_counts, key=type_counts.__getitem__, default="unknown")
        return self._pass(f"Gene IDs appear consistent ({dominant} namespace).")


class DuplicateGeneRule(BaseRule):
    rule_id = "GEN-002"
    category = "gene"
    severity = "ERROR"
    description = "Gene IDs in count matrix must be unique."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        idx = [str(i) for i in context.count_matrix.index]
        seen: dict = {}
        for g in idx:
            seen[g] = seen.get(g, 0) + 1
        dupes = [g for g, cnt in seen.items() if cnt > 1]
        if dupes:
            return self._fail(
                f"{len(dupes)} duplicate gene ID(s) found.",
                affected_items=dupes[:10],
                suggestion=(
                    "Aggregate duplicate genes (e.g. sum counts) or remove duplicates "
                    "before running DE analysis."
                ),
            )
        return self._pass("All gene IDs are unique.")


class VersionSuffixRule(BaseRule):
    rule_id = "GEN-003"
    category = "gene"
    severity = "WARNING"
    description = "Ensembl gene IDs should not include version suffixes (e.g. ENSG00000141510.14)."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        versioned = [
            str(i) for i in context.count_matrix.index
            if _ENSEMBL_VERSION_RE.match(str(i))
        ]
        if versioned:
            return self._fail(
                f"{len(versioned)} Ensembl gene ID(s) have version suffixes.",
                affected_items=versioned[:10],
                suggestion=(
                    "Strip version suffixes (e.g. replace 'ENSG00000141510.14' with "
                    "'ENSG00000141510') for compatibility with annotation databases."
                ),
            )
        return self._pass("No Ensembl version suffixes detected.")


class NonBiologicalIDRule(BaseRule):
    rule_id = "GEN-004"
    category = "gene"
    severity = "WARNING"
    description = "Gene IDs should look biologically meaningful, not random strings."

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")
        gene_ids = [str(i) for i in context.count_matrix.index[:500]]
        unknown = [g for g in gene_ids if _classify_id(g) == "unknown"]
        if len(unknown) > len(gene_ids) * 0.1:
            return self._fail(
                f"{len(unknown)} gene ID(s) do not match any known biological format.",
                affected_items=unknown[:10],
                suggestion=(
                    "Verify that the first column contains gene IDs "
                    "(Ensembl, Entrez, or HGNC gene symbols)."
                ),
            )
        return self._pass("Gene IDs appear biologically meaningful.")
