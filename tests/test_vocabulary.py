from __future__ import annotations

from pathlib import Path

import pytest

import pfpred.vocabulary as vocabulary_module
from pfpred.records import Annotation
from pfpred.vocabulary import (
    VocabularyError,
    _validate_distinct_paths,
    build_vocabulary,
    build_vocabulary_files,
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
    assert coverage["F"]["splits"]["train"]["annotation_pairs"] == 8
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
    annotations = [row(f"T{i}", "GO:1") for i in range(1, 5)]
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


@pytest.mark.parametrize(
    "collision", [
        "annotations_split", "annotations_fasta", "annotations_vocabulary", "annotations_coverage",
        "split_fasta", "split_vocabulary", "split_coverage",
        "fasta_vocabulary", "fasta_coverage", "vocabulary_coverage",
    ]
)
def test_all_input_and_output_path_collisions_are_rejected(
    tmp_path: Path, collision: str,
) -> None:
    paths = {
        "annotations": tmp_path / "annotations.tsv.gz",
        "split": tmp_path / "split.json",
        "fasta": tmp_path / "human.fasta.gz",
        "vocabulary": tmp_path / "vocab.json",
        "coverage": tmp_path / "coverage.json",
    }
    first, second = collision.split("_", maxsplit=1)
    paths[second] = paths[first]
    with pytest.raises(VocabularyError, match="path collision"):
        _validate_distinct_paths(
            paths["annotations"], paths["split"], paths["vocabulary"], paths["coverage"],
            paths["fasta"],
        )


def test_build_vocabulary_files_rejects_output_overwriting_fasta(
    tmp_path: Path,
) -> None:
    annotations_path = tmp_path / "annotations.tsv.gz"
    split_path = tmp_path / "human_split.json"
    fasta_path = tmp_path / "human.fasta.gz"
    coverage_path = tmp_path / "coverage.json"
    fasta_path.write_bytes(b"original raw FASTA")
    for path in (annotations_path, split_path):
        path.write_bytes(b"fixture")

    with pytest.raises(VocabularyError, match="path collision"):
        build_vocabulary_files(
            annotations_path, split_path, fasta_path, coverage_path,
            fasta_path=fasta_path,
        )

    assert fasta_path.read_bytes() == b"original raw FASTA"


def test_unknown_annotation_aspect_raises_vocabulary_error() -> None:
    with pytest.raises(VocabularyError, match="unknown aspect"):
        build_vocabulary([row("P1", "GO:1", aspect="X")], {"P1": "train"}, min_count=1)



def test_build_vocabulary_files_writes_outputs_and_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    annotations_path = tmp_path / "annotations.tsv.gz"
    split_path = tmp_path / "human_split.json"
    fasta_path = tmp_path / "human.fasta.gz"
    vocabulary_path = tmp_path / "human_vocabulary.json"
    coverage_path = tmp_path / "human_coverage.json"
    for path in (annotations_path, split_path, fasta_path):
        path.write_bytes(b"fixture")

    assignments = {f"P{i}": "train" for i in range(5)}
    split_payload = {
        "track": "A",
        "dataset": "fixture",
        "input_fasta": fasta_path.name,
        "input_sha256": "a" * 64,
        "platform": {"system": "linux"},
        "mmseqs": {"version": "fixture"},
        "assignments": {
            protein: {"split": split} for protein, split in assignments.items()
        },
    }
    rows = [
        row(f"P{i}", term, aspect)
        for aspect, term in (("F", "GO:F"), ("P", "GO:P"), ("C", "GO:C"))
        for i in range(5)
    ]
    monkeypatch.setattr(vocabulary_module, "load_split", lambda *_args, **_kwargs: split_payload)
    monkeypatch.setattr(vocabulary_module, "read_annotations", lambda _path: iter(rows))
    monkeypatch.setattr(vocabulary_module, "sha256_file", lambda _path: "b" * 64)

    vocab_payload, coverage_payload = build_vocabulary_files(
        annotations_path, split_path, vocabulary_path, coverage_path,
        fasta_path=fasta_path, min_count=5,
    )

    assert set(vocab_payload["aspects"]) == {"F", "P", "C"}
    assert vocab_payload["provenance"]["input_fasta_sha256"] == "a" * 64
    assert coverage_payload["unassigned_curated_annotation_proteins"] == 0
    assert vocabulary_path.exists() and coverage_path.exists()


def test_build_vocabulary_files_rejects_empty_aspect(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    annotations_path = tmp_path / "annotations.tsv.gz"
    split_path = tmp_path / "human_split.json"
    fasta_path = tmp_path / "human.fasta.gz"
    vocabulary_path = tmp_path / "vocabulary.json"
    coverage_path = tmp_path / "coverage.json"
    for path in (annotations_path, split_path, fasta_path):
        path.write_bytes(b"fixture")
    split_payload = {
        "track": "A",
        "assignments": {
            f"P{i}": {"split": "train"} for i in range(5)
        },
    }
    rows = [row(f"P{i}", "GO:F") for i in range(5)]
    monkeypatch.setattr(vocabulary_module, "load_split", lambda *_args, **_kwargs: split_payload)
    monkeypatch.setattr(vocabulary_module, "read_annotations", lambda _path: iter(rows))
    with pytest.raises(VocabularyError, match="vocabulary is empty"):
        build_vocabulary_files(
            annotations_path, split_path, vocabulary_path, coverage_path,
            fasta_path=fasta_path, min_count=5,
        )
    assert not vocabulary_path.exists()
    assert not coverage_path.exists()



def test_vocab_cli_passes_layout_fasta_to_builder(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str],
) -> None:
    from pfpred import cli

    captured: dict[str, object] = {}

    def fake_builder(
        annotation_path: Path,
        split_path: Path,
        vocabulary_path: Path,
        coverage_path: Path,
        *,
        fasta_path: Path,
        min_count: int,
    ) -> tuple[dict[str, object], dict[str, object]]:
        captured.update({
            "annotation": annotation_path,
            "split": split_path,
            "vocabulary": vocabulary_path,
            "coverage": coverage_path,
            "fasta": fasta_path,
            "min_count": min_count,
        })
        return (
            {"aspects": {aspect: {"term_count": 1} for aspect in ("F", "P", "C")}},
            {"aspects": {}},
        )

    monkeypatch.setattr(cli, "build_vocabulary_files", fake_builder)
    assert cli.main(["vocab", "--root", str(tmp_path)]) == 0
    assert captured["fasta"] == tmp_path / "data/raw/uniprot/human_proteins_2026_04.fasta.gz"
    assert captured["min_count"] == 5
    assert "Wrote vocabulary" in capsys.readouterr().out


def test_numpy_integer_min_count_is_supported() -> None:
    import numpy as np

    annotations = [row(f"T{i}", "GO:1") for i in range(5)]
    vocabularies, _, _ = build_vocabulary(
        annotations, {f"T{i}": "train" for i in range(5)}, min_count=np.int64(5)
    )
    assert vocabularies["F"] == {"GO:1": 5}
