"""Gene ID validation rules."""
from __future__ import annotations

import re

import pandas as pd

from app.models.context import ValidationContext
from app.models.rule_result import RuleResult
from app.rules.base import BaseRule

_ENSEMBL_RE = re.compile(r'^ENS[A-Z]*G\d{11}(\.\d+)?$', re.IGNORECASE)
_ENSEMBL_BASE_RE = re.compile(r'^ENS[A-Z]*G\d{11}$', re.IGNORECASE)
_ENSEMBL_VERSION_RE = re.compile(r'^ENS[A-Z]*G\d{11}\.\d+$', re.IGNORECASE)
# Multi-species Ensembl ID: ENSG (human), ENSMUSG (mouse), ENSRNOG (rat), ENSDARG (zebrafish), etc.
_ENSEMBL_ANY_RE = re.compile(r'^ENS[A-Z]*G\d', re.IGNORECASE)
_ENTREZ_RE = re.compile(r'^\d+$')

# Ensembl species prefixes used by GEN-005
_SPECIES_PREFIXES: dict[str, tuple[str, re.Pattern]] = {
    "ENSG":     ("human",     re.compile(r'^ENSG\d',     re.IGNORECASE)),
    "ENSMUSG":  ("mouse",     re.compile(r'^ENSMUSG\d',  re.IGNORECASE)),
    "ENSRNOG":  ("rat",       re.compile(r'^ENSRNOG\d',  re.IGNORECASE)),
    "ENSDARG":  ("zebrafish", re.compile(r'^ENSDARG\d',  re.IGNORECASE)),
    "ENSGALG":  ("chicken",   re.compile(r'^ENSGALG\d',  re.IGNORECASE)),
}
_SYMBOL_RE = re.compile(r'^[A-Z][A-Z0-9\-]{1,}$')


def _classify_id(gene_id: str) -> str:
    g = gene_id.strip()
    if _ENSEMBL_BASE_RE.match(g):
        return "ensembl"
    if _ENSEMBL_VERSION_RE.match(g):
        return "ensembl_versioned"
    # Multi-species Ensembl IDs (ENSMUSG, ENSRNOG, ENSDARG, etc.)
    if _ENSEMBL_ANY_RE.match(g):
        return "ensembl"
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


class OrganismDetectionRule(BaseRule):
    """GEN-005 — Detect organism from Ensembl gene ID prefixes and store in context flags.

    Biological importance: Several downstream rules (BIO-005 MT thresholds, BIO-006
    housekeeping gene sets) require knowledge of the organism.  Without explicit
    detection, those rules must either use wrong thresholds or SKIP entirely.
    Additionally, mixed-organism gene IDs (e.g. mouse IDs in a human experiment)
    indicate a critical data error that no other rule would catch.

    The Ensembl species prefix is a stable, published identifier — it is unambiguous
    and does not require any network access or external database lookup.
    """

    rule_id = "GEN-005"
    category = "gene"
    severity = "INFO"
    description = (
        "Detect organism from Ensembl gene ID prefixes and store in validation context. "
        "Supports human (ENSG), mouse (ENSMUSG), rat (ENSRNOG), zebrafish (ENSDARG), "
        "and chicken (ENSGALG). Mixed-organism IDs are flagged as ERROR."
    )

    # Fraction of IDs that must match a single species to call it confidently
    _DOMINANT_THRESHOLD = 0.80
    # Minimum fraction for a species to be considered "present" (cross-contamination check)
    _MINOR_THRESHOLD = 0.10

    def run(self, context: ValidationContext) -> RuleResult:
        if context.count_matrix is None:
            return self._skip("Count matrix not available.")

        gene_ids = [str(g) for g in context.count_matrix.index]
        n = len(gene_ids)
        if n == 0:
            return self._skip("Count matrix has no rows.")

        # Count how many IDs match each species prefix
        species_counts: dict[str, int] = {}
        for label, (species_name, pattern) in _SPECIES_PREFIXES.items():
            cnt = sum(1 for g in gene_ids if pattern.match(g))
            if cnt > 0:
                species_counts[species_name] = cnt

        if not species_counts:
            return self._skip(
                "No Ensembl species prefix recognised; organism detection skipped. "
                "Gene IDs may be Entrez IDs or gene symbols."
            )

        total_ensembl = sum(species_counts.values())
        fracs = {sp: cnt / n for sp, cnt in species_counts.items()}

        # Check for mixed-organism IDs (two species each > 10%)
        present_species = [sp for sp, frac in fracs.items() if frac > self._MINOR_THRESHOLD]
        if len(present_species) > 1:
            detail_str = ", ".join(
                f"{sp}: {fracs[sp]*100:.1f}%" for sp in sorted(present_species)
            )
            return self._fail(
                f"Mixed-organism Ensembl gene IDs detected: {detail_str}. "
                "Gene IDs from multiple species in a single matrix indicate a data error.",
                affected_items=[f"{sp}: {species_counts[sp]} IDs" for sp in sorted(present_species)],
                suggestion=(
                    "Ensure all gene IDs belong to a single organism. Remove or separate "
                    "gene rows from different species."
                ),
                details={"species_fractions": fracs},
            )

        # Single dominant organism
        dominant_species, dominant_frac = max(fracs.items(), key=lambda x: x[1])
        if dominant_frac >= self._DOMINANT_THRESHOLD:
            context.flags["organism"] = dominant_species
            return RuleResult(
                rule_id=self.rule_id,
                category=self.category,
                severity=self.severity,
                status="PASS",
                message=(
                    f"Organism detected: {dominant_species!r} "
                    f"({dominant_frac*100:.1f}% of gene IDs match Ensembl prefix)."
                ),
                details={"organism": dominant_species, "species_fractions": fracs},
            )

        # Some Ensembl IDs but below confidence threshold
        context.flags["organism"] = dominant_species
        return self._skip(
            f"Organism detection inconclusive: best match is {dominant_species!r} "
            f"({dominant_frac*100:.1f}% < {self._DOMINANT_THRESHOLD*100:.0f}% threshold). "
            "Using best guess for downstream rules."
        )
