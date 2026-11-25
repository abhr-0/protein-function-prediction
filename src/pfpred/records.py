"""Records shared by every adapter."""

from __future__ import annotations

from dataclasses import dataclass

VALID_ASPECTS = frozenset({"F", "P", "C"})


class DataFormatError(ValueError):
    """Raised when a raw data file does not match its documented format."""


@dataclass(frozen=True, slots=True)
class Protein:
    """A protein sequence keyed by its primary accession."""

    accession: str
    sequence: str


@dataclass(frozen=True, slots=True)
class Annotation:
    """One protein-term annotation.

    ``evidence`` is a UniProt evidence code (``IDA``, ``IEA``, ...). CAFA6
    labels carry no evidence code, so their ``evidence`` is empty.
    """

    entry: str
    term: str
    aspect: str
    evidence: str = ""


def validate_aspect(aspect: str, *, source: str) -> str:
    """Return ``aspect`` if it is one of F/P/C, else raise :class:`DataFormatError`."""
    if aspect not in VALID_ASPECTS:
        raise DataFormatError(f"{source}: unknown aspect {aspect!r}")
    return aspect
