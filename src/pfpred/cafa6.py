"""CAFA6 adapters (Track B).

CAFA6 labels carry no evidence code, so every :class:`Annotation` produced here
has ``evidence=""``.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from pathlib import Path

from pfpred.io import expect_header, read_tsv
from pfpred.records import Annotation, DataFormatError, Protein, validate_aspect
from pfpred.uniprot import read_fasta

TERMS_HEADER = ("EntryID", "term", "aspect")


def _read_term_rows(path: Path, header: tuple[str, ...]) -> Iterator[Annotation]:
    rows = read_tsv(path)
    try:
        first = next(rows)
    except StopIteration:
        raise DataFormatError(f"{path}: empty file") from None
    expect_header(first, header, source=str(path))
    for row in rows:
        if len(row) != len(header):
            raise DataFormatError(f"{path}: expected {len(header)} fields, got {row!r}")
        entry, term, aspect = (field.strip() for field in row)
        if not entry or not term:
            raise DataFormatError(f"{path}: empty entry or term in {row!r}")
        validate_aspect(aspect, source=str(path))
        yield Annotation(entry=entry, term=term, aspect=aspect)


def read_terms(path: Path) -> Iterator[Annotation]:
    """Read CAFA6 ``train_terms.tsv`` (t0 labels)."""
    return _read_term_rows(path, TERMS_HEADER)


def read_groundtruth(path: Path) -> Iterator[Annotation]:
    """Read one ``May2026_groundtruth_*.tsv`` file."""
    return _read_term_rows(path, TERMS_HEADER)


def read_taxonomy(path: Path) -> dict[str, int]:
    """Read the headerless ``train_taxonomy.tsv`` as ``protein -> taxid``."""
    taxonomy: dict[str, int] = {}
    for row in read_tsv(path):
        if len(row) != 2:
            raise DataFormatError(f"{path}: expected 2 fields, got {row!r}")
        entry, raw_taxid = (field.strip() for field in row)
        if not entry or not raw_taxid:
            raise DataFormatError(f"{path}: empty entry or taxid in {row!r}")
        try:
            taxid = int(raw_taxid)
        except ValueError:
            raise DataFormatError(f"{path}: non-integer taxid {raw_taxid!r}") from None
        if taxid <= 0:
            raise DataFormatError(f"{path}: non-positive taxid {taxid}")
        previous = taxonomy.setdefault(entry, taxid)
        if previous != taxid:
            raise DataFormatError(f"{path}: conflicting taxids for {entry}: {previous}, {taxid}")
    return taxonomy


def read_sequences(path: Path) -> Iterator[Protein]:
    """Read CAFA6 sequences, accepting both header styles."""
    return read_fasta(path)


def groundtruth_proteins(path: Path) -> frozenset[str]:
    return frozenset(row.entry for row in read_groundtruth(path))


def groundtruth_union(paths: Iterable[Path]) -> frozenset[str]:
    """Union the proteins of the given ground-truth files."""
    proteins: set[str] = set()
    for path in paths:
        proteins |= groundtruth_proteins(path)
    return frozenset(proteins)
