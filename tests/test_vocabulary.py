from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pytest

from pfpred.records import Annotation
from pfpred.vocabulary import (
    VocabularyError,
    _validate_distinct_paths,
    build_vocabulary,
)


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
        "T5": "train", "T6": "train", "V1": "validation", "X1": "train",
    }

    vocabularies, coverage, unassigned = build_vocabulary(
        annotations, assignments, min_count=5
    )

    assert unassigned == 0
    assert vocabularies["F"] == {"GO:1": 5}
    assert vocabularies["P"] == {}
    assert vocabularies["C"] == {}
    assert coverage["F"]["terms_kept"] == 1
    assert coverage["F"]["splits"]["validation"]["annotation_pairs_covered"] == 1
    assert coverage["F"]["splits"]["validation"]["proteins_dropped_by_min_count"] == 0
    assert coverage["F"]["splits"]["train"]["annotation_pairs"] == 7
    assert coverage["F"]["splits"]["train"]["proteins_in_split"] == 7
    assert coverage["F"]["splits"]["train"]["proteins_without_aspect_label"] == 1
    assert coverage["F"]["splits"]["train"]["proteins_dropped_by_min_count"] == 1


def test_holdout_annotations_cannot_enter_vocabulary() -> None:
    annotations = [row(f"T{i}", "GO:train") for i in range(5)]
    assignments = {f"T{i}": "train" for i in range(5)}
    for i in range(5):
        annotations.append(row(f"V{i}", "GO:leak_val"))
        assignments[f"V{i}"] = "validation"
    for i in range(5):
        annotations.append(row(f"E{i}", "GO:leak_test"))
        assignments[f"E{i}"] = "test"

    vocabularies, _, _ = build_vocabulary(annotations, assignments, min_count=5)
    assert vocabularies["F"] == {"GO:train": 5}


def test_frequency_boundary_and_duplicate_rows_are_exact() -> None:
    annotations = [row(f"T{i}", "GO:kept") for i in range(5)]
    annotations.extend(row("T0", "GO:kept", evidence="IMP") for _ in range(3))
    annotations.extend(row(f"T{i}", "GO:dropped") for i in range(4))
    assignments = {f"T{i}": "train" for i in range(5)}

    vocabularies, _, _ = build_vocabulary(annotations, assignments, min_count=5)
    assert vocabularies["F"] == {"GO:kept": 5}


def test_selected_terms_have_a_training_positive() -> None:
    annotations = [row(f"T{i}", "GO:1") for i in range(5)]
    vocabularies, coverage, _ = build_vocabulary(
        annotations, {f"T{i}": "train" for i in range(5)}, min_count=5
    )
    assert vocabularies["F"]["GO:1"] == 5
    assert coverage["F"]["train_positive_proteins"] == 5


def test_curated_policy_is_applied_before_counting() -> None:
    annotations = [row(f"T{i}", "GO:1") for i in range(4)]
    annotations.append(row("T5", "GO:1", evidence="IEA"))
    vocabularies, _, _ = build_vocabulary(
        annotations, {f"T{i}": "train" for i in range(1, 6)}, min_count=5
    )
    assert "GO:1" not in vocabularies["F"]


def test_unassigned_curated_protein_is_a_hard_error() -> None:
    annotations = [row("OUTSIDE", "GO:1")]
    with pytest.raises(VocabularyError, match="absent from the split"):
        build_vocabulary(annotations, {}, min_count=1)


@pytest.mark.parametrize("min_count", [0, -1, True, 1.5])
def test_min_count_must_be_positive_integer(min_count: int) -> None:
    with pytest.raises(VocabularyError, match="min_count"):
        build_vocabulary([], {}, min_count=min_count)


def test_path_collisions_are_rejected(tmp_path: Path) -> None:
    annotations = tmp_path / "annotations.tsv.gz"
    split = tmp_path / "split.json"
    output = tmp_path / "vocab.json"
    with pytest.raises(VocabularyError, match="path collision"):
        _validate_distinct_paths(annotations, split, output, output)
