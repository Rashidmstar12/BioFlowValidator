"""FileParser: reads and normalises uploaded files into ValidationContext."""
from __future__ import annotations

import io
from typing import Optional, Tuple

import chardet
import pandas as pd

from app.models.context import ValidationContext


def _detect_encoding(data: bytes) -> str:
    result = chardet.detect(data)
    enc = result.get("encoding") or "utf-8"
    # Normalise common aliases
    enc = enc.lower().replace("-", "")
    if enc in ("utf8", "utf8bom", "utf8sig"):
        return "utf-8-sig"
    return result.get("encoding") or "utf-8"


def _detect_delimiter(text: str) -> str:
    """Return the most likely delimiter from the first two lines."""
    lines = [l for l in text.splitlines() if l.strip()][:2]
    sample = "\n".join(lines)
    tab_count = sample.count("\t")
    comma_count = sample.count(",")
    return "\t" if tab_count >= comma_count else ","


def _read_tabular(
    data: bytes, filename: str
) -> Tuple[Optional[pd.DataFrame], str, str]:
    """Attempt to parse bytes as a tabular file. Returns (df, delimiter, encoding)."""
    encoding = _detect_encoding(data)
    try:
        text = data.decode(encoding, errors="replace")
    except Exception:
        text = data.decode("utf-8", errors="replace")

    delimiter = _detect_delimiter(text)

    try:
        df = pd.read_csv(
            io.StringIO(text),
            sep=delimiter,
            index_col=0,
            engine="python",
        )
        return df, delimiter, encoding
    except Exception:
        return None, delimiter, encoding


def parse_files(
    count_bytes: bytes,
    count_filename: str,
    metadata_bytes: Optional[bytes] = None,
    metadata_filename: str = "",
) -> ValidationContext:
    """Parse uploaded files and return a populated ValidationContext."""
    ctx = ValidationContext(
        count_matrix_bytes=count_bytes,
        count_filename=count_filename,
        metadata_bytes=metadata_bytes or b"",
        metadata_filename=metadata_filename,
    )

    count_df, count_delim, _ = _read_tabular(count_bytes, count_filename)
    ctx.count_matrix = count_df
    ctx.count_delimiter = count_delim

    if metadata_bytes:
        meta_df, meta_delim, _ = _read_tabular(metadata_bytes, metadata_filename)
        ctx.metadata = meta_df
        ctx.metadata_delimiter = meta_delim

    return ctx
