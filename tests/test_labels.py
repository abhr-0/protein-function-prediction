from __future__ import annotations

import pytest

from pfpred.labels import (
    CURATED,
    EvidencePolicy,
    filter_annotations,
    group_by_aspect,
    group_by_protein,
    summarize,
)
from pfpred.records import Annotation, DataFormatError


def annotation(entry: str, term: str, aspect: str = "F", evidence: str = "IDA") -> Annotation:
    return Annotation(entry=entry, term=term, aspect=aspect, evidence=evidence)


@pytest.mark.parametrize("evidence", ["IEA", "ND", "NR"])
def test_curated_drops_excluded_evidence(evidence: str) -> None:
    assert not CURATED.keeps(evidence)


@pytest.mark.parametrize("evidence", ["IDA", "IPI", "IBA", "ISS", "TAS", "EXP", "ISM", ""])
def test_curated_keeps_other_evidence(evidence: str) -> None:
    assert CURATED.keeps(evidence)


def test_custom_policy_can_exclude_more() -> None:
    strict = EvidencePolicy("strict", frozenset({"IEA"}))
    assert strict.keeps("ND")
    assert not strict.keeps("IEA")


def test_filter_annotations() -> None:
    rows = [
        annotation("P1", "GO:1", evidence="IDA"),
        annotation("P1", "GO:2", evidence="IEA"),
        annotation("P2", "GO:3", evidence="ND"),
        annotation("P2", "GO:4", evidence="TAS"),
    ]
    kept = list(filter_annotations(rows))
    assert [(a.entry, a.term) for a in kept] == [("P1", "GO:1"), ("P2", "GO:4")]


def test_group_by_protein_folds_duplicates() -> None:
    rows = [
        annotation("P1", "GO:1", evidence="IDA"),
        annotation("P1", "GO:1", evidence="IMP"),
        annotation("P1", "GO:2"),
        annotation("P2", "GO:3"),
    ]
    grouped = group_by_protein(rows)
    assert grouped == {"P1": frozenset({"GO:1", "GO:2"}), "P2": frozenset({"GO:3"})}


def test_group_by_aspect_separates_aspects() -> None:
    rows = [
        annotation("P1", "GO:1", aspect="F"),
        annotation("P1", "GO:2", aspect="P"),
        annotation("P2", "GO:3", aspect="C"),
    ]
    grouped = group_by_aspect(rows)
    assert grouped == {
        "F": {"P1": frozenset({"GO:1"})},
        "P": {"P1": frozenset({"GO:2"})},
        "C": {"P2": frozenset({"GO:3"})},
    }


def test_group_by_protein_rejects_bad_aspect() -> None:
    with pytest.raises(DataFormatError, match="unknown aspect"):
        group_by_protein([annotation("P1", "GO:1", aspect="X")])


def test_summarize() -> None:
    rows = [
        annotation("P1", "GO:1", aspect="F"),
        annotation("P1", "GO:2", aspect="F"),
        annotation("P2", "GO:3", aspect="P"),
    ]
    stats = summarize(rows)
    assert stats.annotations == 3
    assert stats.proteins == 2
    assert stats.terms == 3
    assert stats.per_aspect == {"F": 2, "P": 1}


def test_summarize_empty() -> None:
    stats = summarize([])
    assert stats.annotations == 0
    assert stats.proteins == 0
    assert stats.terms == 0
    assert stats.per_aspect == {}
