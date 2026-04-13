"""Rule registry — central list of all validation rules in execution order."""
from __future__ import annotations

from app.rules.base import BaseRule
from app.rules.format.rules import (
    EncodingRule,
    DelimiterRule,
    HeaderRule,
    DuplicateColumnRule,
    NonNumericRule,
    NegativeCountRule,
    WhitespaceNameRule,
    MatrixOrientationRule,
)
from app.rules.sample.rules import (
    SampleMatchRule,
    SampleOrderRule,
    DuplicateSampleRule,
    ReplicateCountRule,
    NearIdenticalSampleRule,
)
from app.rules.gene.rules import (
    GeneIDFormatRule,
    DuplicateGeneRule,
    VersionSuffixRule,
    NonBiologicalIDRule,
    OrganismDetectionRule,
)
from app.rules.normalization.rules import (
    NonIntegerCountRule,
    LibrarySizeRule,
    ZeroLibraryRule,
    AllZeroGeneRule,
    LowCountDominanceRule,
    DuplicateRowRule,
)
from app.rules.biology.rules import (
    SingleConditionRule,
    ConditionLabelSanityRule,
    MetadataCardinalityRule,
    HighCountGeneRule,
    MitochondrialFractionRule,
    HousekeepingGeneRule,
    BatchConfoundingRule,
    ERCCSpikeInRule,
)

# Rules are executed in this order.
# Format rules run first (gate for downstream checks).
# GEN-005 (OrganismDetectionRule) runs early in the gene block so that downstream
# biology rules (BIO-005, BIO-006) can read context.flags["organism"].
RULE_CLASSES: list[type[BaseRule]] = [
    # Format
    EncodingRule,
    DelimiterRule,
    MatrixOrientationRule,   # FMT-008: check orientation before parsing content
    HeaderRule,
    DuplicateColumnRule,
    NonNumericRule,
    NegativeCountRule,
    WhitespaceNameRule,
    # Sample
    DuplicateSampleRule,
    SampleMatchRule,
    SampleOrderRule,
    ReplicateCountRule,
    NearIdenticalSampleRule,
    # Gene
    OrganismDetectionRule,   # GEN-005: sets context.flags["organism"] for BIO-005/006
    GeneIDFormatRule,
    DuplicateGeneRule,
    VersionSuffixRule,
    NonBiologicalIDRule,
    # Normalization
    NonIntegerCountRule,
    LibrarySizeRule,
    ZeroLibraryRule,
    AllZeroGeneRule,
    LowCountDominanceRule,
    DuplicateRowRule,
    # Biology
    SingleConditionRule,
    ConditionLabelSanityRule,
    MetadataCardinalityRule,
    HighCountGeneRule,
    MitochondrialFractionRule,
    HousekeepingGeneRule,
    BatchConfoundingRule,    # BIO-007
    ERCCSpikeInRule,         # BIO-008
]


def get_rules() -> list[BaseRule]:
    """Return instantiated rule objects in execution order."""
    return [cls() for cls in RULE_CLASSES]
