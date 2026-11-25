"""Evidence policy and label grouping.

Training uses raw annotations; ancestors are only added at evaluation (Track B,
via ``cafaeval -prop fill``). This module therefore does not read any ontology.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from pfpred.records import Annotation, validate_aspect

EXCLUDED_EVIDENCE = frozenset({"IEA", "ND", "NR"})


@dataclass(frozen=True, slots=True)
class EvidencePolicy:
    """Keep annotations whose evidence code is not excluded."""

    name: str
    excluded: frozenset[str] = EXCLUDED_EVIDENCE

    def keeps(self, evidence: str) -> bool:
        return evidence not in self.excluded


CURATED = EvidencePolicy("curated")
"""The project's Track A policy: drop IEA, ND, and NR."""


def filter_annotations(
    annotations: Iterable[Annotation],
    policy: EvidencePolicy = CURATED,
) -> Iterator[Annotation]:
    for annotation in annotations:
        if policy.keeps(annotation.evidence):
            yield annotation


def group_by_protein(annotations: Iterable[Annotation]) -> dict[str, frozenset[str]]:
    """Map each protein to its set of terms, validating aspects on the way."""
    grouped: dict[str, set[str]] = defaultdict(set)
    for annotation in annotations:
        validate_aspect(annotation.aspect, source="annotations")
        grouped[annotation.entry].add(annotation.term)
    return {entry: frozenset(terms) for entry, terms in grouped.items()}


def group_by_aspect(
    annotations: Iterable[Annotation],
) -> dict[str, dict[str, frozenset[str]]]:
    """Split annotations by aspect, then group each aspect by protein."""
    by_aspect: dict[str, list[Annotation]] = defaultdict(list)
    for annotation in annotations:
        validate_aspect(annotation.aspect, source="annotations")
        by_aspect[annotation.aspect].append(annotation)
    return {aspect: group_by_protein(rows) for aspect, rows in by_aspect.items()}


@dataclass(frozen=True, slots=True)
class LabelStats:
    annotations: int
    proteins: int
    terms: int
    per_aspect: dict[str, int]


def summarize(annotations: Iterable[Annotation]) -> LabelStats:
    """Count rows, proteins, terms, and rows per aspect."""
    rows = 0
    proteins: set[str] = set()
    terms: set[str] = set()
    per_aspect: dict[str, int] = defaultdict(int)
    for annotation in annotations:
        validate_aspect(annotation.aspect, source="annotations")
        rows += 1
        proteins.add(annotation.entry)
        terms.add(annotation.term)
        per_aspect[annotation.aspect] += 1
    return LabelStats(
        annotations=rows,
        proteins=len(proteins),
        terms=len(terms),
        per_aspect=dict(per_aspect),
    )
