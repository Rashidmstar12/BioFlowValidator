"""ValidationContext holds all parsed data passed to rules."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import pandas as pd


@dataclass
class ValidationContext:
    """Container for all parsed input data made available to every rule."""

    # Raw bytes of uploaded files (for encoding checks)
    count_matrix_bytes: bytes = b""
    metadata_bytes: bytes = b""

    # Parsed DataFrames (genes × samples orientation; index = gene IDs)
    count_matrix: Optional[pd.DataFrame] = None
    # Metadata: index = sample IDs, columns = phenotype variables
    metadata: Optional[pd.DataFrame] = None

    # Detected delimiter for each file
    count_delimiter: str = ""
    metadata_delimiter: str = ""

    # Filename strings (for display only)
    count_filename: str = ""
    metadata_filename: str = ""

    # Extra flags set by earlier rules so later rules can act on them
    flags: dict = field(default_factory=dict)
