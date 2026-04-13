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
)

# Rules are executed in this order.
# Format rules run first (gate for downstream checks).
RULE_CLASSES: list[type[BaseRule]] = [
    # Format
    EncodingRule,
    DelimiterRule,
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
]


def get_rules() -> list[BaseRule]:
    """Return instantiated rule objects in execution order."""
    return [cls() for cls in RULE_CLASSES]
