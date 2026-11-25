"""Small I/O helpers shared by the adapters."""

from __future__ import annotations

import csv
import gzip
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import TextIO

from pfpred.records import DataFormatError


def open_text(path: Path) -> TextIO:
    """Open a UTF-8 text file, transparently decompressing ``.gz``."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, encoding="utf-8", newline="")


def read_tsv(path: Path, *, delimiter: str = "\t") -> Iterator[list[str]]:
    """Yield non-empty rows of a TSV file as lists of fields.

    A header row, if present, is yielded like any other row; callers validate it
    with :func:`expect_header`.
    """
    with open_text(path) as handle:
        for row in csv.reader(handle, delimiter=delimiter):
            if row:
                yield row


def expect_header(row: Sequence[str], expected: Sequence[str], *, source: str) -> None:
    """Raise :class:`DataFormatError` unless ``row`` matches ``expected``."""
    if list(row) != list(expected):
        raise DataFormatError(f"{source}: expected header {list(expected)}, got {list(row)}")
