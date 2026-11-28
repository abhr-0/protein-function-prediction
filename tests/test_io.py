from __future__ import annotations

from pathlib import Path

import pytest

from pfpred.io import expect_header, open_text, read_tsv
from pfpred.records import DataFormatError


def test_read_tsv_plain(tmp_path: Path) -> None:
    path = tmp_path / "tiny.tsv"
    path.write_text("a\tb\n\n1\t2\n", encoding="utf-8")
    assert list(read_tsv(path)) == [["a", "b"], ["1", "2"]]


def test_read_tsv_gzip(tmp_path: Path, gzipped) -> None:
    source = tmp_path / "tiny.tsv"
    source.write_text("a\tb\n1\t2\n", encoding="utf-8")
    packed = gzipped(source)
    assert list(read_tsv(packed)) == [["a", "b"], ["1", "2"]]


def test_open_text_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        with open_text(tmp_path / "absent.tsv") as handle:
            handle.read()


def test_expect_header_accepts_match() -> None:
    expect_header(["Entry", "GO"], ["Entry", "GO"], source="tiny")


def test_expect_header_rejects_mismatch() -> None:
    with pytest.raises(DataFormatError, match="tiny: expected header"):
        expect_header(["Entry", "Term"], ["Entry", "GO"], source="tiny")


def test_write_json_atomic_is_deterministic_and_cleans_temporary_files(
    tmp_path: Path,
) -> None:
    from pfpred.io import write_json_atomic

    output = tmp_path / "nested" / "result.json"
    payload = {"z": 2, "a": {"value": 1}}
    write_json_atomic(output, payload)
    first = output.read_bytes()
    write_json_atomic(output, payload)
    assert output.read_bytes() == first
    assert first.endswith(b"\\n")
    assert list(output.parent.iterdir()) == [output]
