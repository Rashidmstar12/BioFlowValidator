"""Unit tests for FileParser."""
from __future__ import annotations

import pytest

from app.engine.parser import parse_files, _detect_delimiter, _detect_encoding


TSV_CONTENT = b"gene_id\tctrl_1\tctrl_2\nENSG00000000001\t10\t20\nENSG00000000002\t5\t15\n"
CSV_CONTENT = b"gene_id,ctrl_1,ctrl_2\nENSG00000000001,10,20\nENSG00000000002,5,15\n"
META_TSV = b"sample_id\tcondition\nctrl_1\tcontrol\nctrl_2\tcontrol\n"


def test_detect_delimiter_tsv():
    assert _detect_delimiter("gene_id\ts1\ts2\nA\t1\t2") == "\t"


def test_detect_delimiter_csv():
    assert _detect_delimiter("gene_id,s1,s2\nA,1,2") == ","


def test_parse_tsv_count_matrix():
    ctx = parse_files(TSV_CONTENT, "counts.tsv")
    assert ctx.count_matrix is not None
    assert list(ctx.count_matrix.columns) == ["ctrl_1", "ctrl_2"]
    assert ctx.count_delimiter == "\t"


def test_parse_csv_count_matrix():
    ctx = parse_files(CSV_CONTENT, "counts.csv")
    assert ctx.count_matrix is not None
    assert list(ctx.count_matrix.columns) == ["ctrl_1", "ctrl_2"]
    assert ctx.count_delimiter == ","


def test_parse_with_metadata():
    ctx = parse_files(TSV_CONTENT, "counts.tsv", META_TSV, "metadata.tsv")
    assert ctx.metadata is not None
    assert list(ctx.metadata.index) == ["ctrl_1", "ctrl_2"]
    assert "condition" in ctx.metadata.columns


def test_parse_empty_bytes():
    ctx = parse_files(b"", "empty.tsv")
    assert ctx.count_matrix is None or ctx.count_matrix.empty


def test_file_hashes_stored():
    ctx = parse_files(TSV_CONTENT, "counts.tsv")
    assert ctx.count_matrix_bytes == TSV_CONTENT
