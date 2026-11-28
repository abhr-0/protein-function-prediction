"""Small I/O helpers shared by the adapters."""

from __future__ import annotations

import csv
import gzip
import json
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import TextIO

from pfpred.records import DataFormatError


def open_text(path: Path) -> TextIO:
    """Open a UTF-8 text file, transparently decompressing ``.gz``."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, encoding="utf-8", newline="")


def read_tsv(path: Path, *, delimiter: str = "\t") -> Iterator[list[str]]:
    """Yield non-empty rows of a TSV file as lists of fields."""
    with open_text(path) as handle:
        for row in csv.reader(handle, delimiter=delimiter):
            if row:
                yield row


def expect_header(row: Sequence[str], expected: Sequence[str], *, source: str) -> None:
    """Raise ``DataFormatError`` unless ``row`` matches ``expected``."""
    if list(row) != list(expected):
        raise DataFormatError(f"{source}: expected header {list(expected)}, got {list(row)}")


def write_json_atomic(path: Path, payload: Mapping[str, object]) -> None:
    """Write sorted, indented JSON through a temporary file in the destination directory."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        temporary.replace(path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
