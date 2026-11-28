from __future__ import annotations

import pytest

from pfpred.records import Annotation
from pfpred.vocabulary import VocabularyError, build_vocabulary


def row(protein: str, term: str, aspect: str = "F", evidence: str = "IDA") -> Annotation:
    return Annotation(protein, term, aspect, evidence)


def test_vocabulary_uses_distinct_training_proteins_and_separates_aspects() -> None:
    annotations = [
        row("T1", "GO:1"), row("T1", "GO:1", evidence="IMP"),
        row("T2", "GO:1"), row("T3", "GO:1"),
        row("T4", "GO:1"), row("T5", "GO:1"),
        row("T1", "GO:2"), row("T2", "GO:2"),
        row("T1", "GO:3", aspect="P"),
        row("V1", "GO:1"), row("X1", "GO:4"),
    ]
    assignments = {
        "T1": "train", "T2": "train", "T3": "train", "T4": "train",
        "T5": "train", "V1": "validation",
    }

    vocabularies, coverage = build_vocabulary(annotations, assignments, min_count=5)

    assert vocabularies["F"] == {"GO:1": 5}
    assert vocabularies["P"] == {}
    assert vocabularies["C"] == {}
    assert coverage["F"]["terms_kept"] == 1
    assert coverage["F"]["splits"]["validation"]["annotation_pairs_covered"] == 1
    assert coverage["F"]["splits"]["validation"]["proteins_dropped_zero_positives"] == 0
    assert coverage["F"]["splits"]["train"]["annotation_pairs"] == 7


def test_holdout_annotations_cannot_change_vocabulary() -> None:
    train = [row(f"T{i}", "GO:train") for i in range(5)]
    assignments = {f"T{i}": "train" for i in range(5)}
    assignments.update({"V1": "validation", "E1": "test"})
    first, _ = build_vocabulary(
        [*train, row("V1", "GO:holdout"), row("E1", "GO:holdout")],
        assignments, min_count=5,
    )
    second, _ = build_vocabulary(
        [*train, row("V1", "GO:other"), row("E1", "GO:other")],
        assignments, min_count=5,
    )
    assert first == second
    assert first["F"] == {"GO:train": 5}


def test_selected_terms_have_a_training_positive() -> None:
    annotations = [row(f"T{i}", "GO:1") for i in range(5)]
    vocabularies, coverage = build_vocabulary(
        annotations, {f"T{i}": "train" for i in range(5)}, min_count=5
    )
    assert vocabularies["F"]["GO:1"] >= 5
    assert coverage["F"]["train_positive_proteins"] == 5


def test_curated_policy_is_applied_before_counting() -> None:
    annotations = [row(f"T{i}", "GO:1") for i in range(4)]
    annotations.append(row("T5", "GO:1", evidence="IEA"))
    vocabularies, _ = build_vocabulary(
        annotations, {f"T{i}": "train" for i in range(1, 6)}, min_count=5
    )
    assert "GO:1" not in vocabularies["F"]


@pytest.mark.parametrize("min_count", [0, -1])
def test_min_count_must_be_positive(min_count: int) -> None:
    with pytest.raises(VocabularyError, match="min_count"):
        build_vocabulary([], {}, min_count=min_count)
