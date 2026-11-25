"""Human UniProt adapters (Track A).

Inputs are the two files derived once from the UniProtKB REST ``format=txt``
stream; see ``docs/data_sources.md``.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from pfpred.io import expect_header, open_text, read_tsv
from pfpred.labels import CURATED, EvidencePolicy, filter_annotations, group_by_protein
from pfpred.records import Annotation, DataFormatError, Protein, validate_aspect

ANNOTATION_HEADER = ("Entry", "GO", "Aspect", "Evidence")


def normalize_accession(header: str) -> str:
    """Extract the primary accession from a FASTA header.

    Accepts a bare accession (``P12345``), a UniProt-style header
    (``sp|P12345|NAME``), and a header with trailing fields
    (``P12345 9606``).
    """
    token = header.strip().split(" ", 1)[0]
    if not token:
        raise DataFormatError("empty FASTA header")
    if "|" in token:
        parts = token.split("|")
        if len(parts) < 2 or not parts[1]:
            raise DataFormatError(f"malformed FASTA header {header!r}")
        return parts[1]
    return token


def read_fasta(path: Path) -> Iterator[Protein]:
    """Stream a FASTA file, normalising accessions and uppercasing sequences."""
    accession: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if accession is not None:
                    yield Protein(accession, "".join(chunks).upper())
                accession = normalize_accession(line[1:])
                chunks = []
            else:
                if accession is None:
                    raise DataFormatError(f"{path}: sequence data before the first header")
                chunks.append("".join(line.split()))
    if accession is not None:
        yield Protein(accession, "".join(chunks).upper())


def read_annotations(path: Path) -> Iterator[Annotation]:
    """Stream ``Entry, GO, Aspect, Evidence`` rows."""
    rows = read_tsv(path)
    try:
        header = next(rows)
    except StopIteration:
        raise DataFormatError(f"{path}: empty annotation file") from None
    expect_header(header, ANNOTATION_HEADER, source=str(path))
    for row in rows:
        if len(row) != len(ANNOTATION_HEADER):
            raise DataFormatError(f"{path}: expected 4 fields, got {row!r}")
        entry, term, aspect, evidence = (field.strip() for field in row)
        if not entry or not term:
            raise DataFormatError(f"{path}: empty entry or term in {row!r}")
        validate_aspect(aspect, source=str(path))
        yield Annotation(entry=entry, term=term, aspect=aspect, evidence=evidence)


def load_labels(
    path: Path,
    *,
    aspect: str,
    policy: EvidencePolicy = CURATED,
) -> dict[str, frozenset[str]]:
    """Return ``protein -> terms`` for one aspect after applying ``policy``."""
    validate_aspect(aspect, source=str(path))
    rows = (
        row
        for row in filter_annotations(read_annotations(path), policy)
        if row.aspect == aspect
    )
    return group_by_protein(rows)
