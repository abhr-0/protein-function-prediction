from __future__ import annotations

from pathlib import Path

import pytest

from pfpred.records import Annotation, DataFormatError
from pfpred.uniprot import load_labels, normalize_accession, read_annotations, read_fasta


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        ("A0A087X1C5", "A0A087X1C5"),
        ("sp|P12345|TEST_HUMAN", "P12345"),
        ("tr|Q99999|Q99999_HUMAN extra words", "Q99999"),
        ("P12345 9606", "P12345"),
    ],
)
def test_normalize_accession(header: str, expected: str) -> None:
    assert normalize_accession(header) == expected


@pytest.mark.parametrize("header", ["", "   ", "sp||NAME"])
def test_normalize_accession_rejects_malformed(header: str) -> None:
    with pytest.raises(DataFormatError, match="FASTA header|empty"):
        normalize_accession(header)


def test_read_fasta(fixtures_dir: Path) -> None:
    proteins = list(read_fasta(fixtures_dir / "human_proteins_mini.fasta"))
    assert [(p.accession, p.sequence) for p in proteins] == [
        ("A0A087X1C5", "MGLEALVPL"),
        ("A0A0C5B5G6", ""),
        ("P12345", "MKV"),
    ]


def test_read_fasta_gzip(fixtures_dir: Path, gzipped) -> None:
    packed = gzipped(fixtures_dir / "human_proteins_mini.fasta")
    assert len(list(read_fasta(packed))) == 3


def test_read_fasta_rejects_sequence_before_header(tmp_path: Path) -> None:
    path = tmp_path / "bad.fasta"
    path.write_text("MKV\n", encoding="utf-8")
    with pytest.raises(DataFormatError, match="before the first header"):
        list(read_fasta(path))


def test_read_annotations(fixtures_dir: Path) -> None:
    rows = list(read_annotations(fixtures_dir / "human_annotations_mini.tsv"))
    assert rows[0] == Annotation(entry="P1", term="GO:0000001", aspect="F", evidence="IDA")
    assert len(rows) == 3


def test_read_annotations_gzip(fixtures_dir: Path, gzipped) -> None:
    packed = gzipped(fixtures_dir / "human_annotations_mini.tsv")
    assert len(list(read_annotations(packed))) == 3


def test_read_annotations_rejects_wrong_header(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("Entry\tGO\tAspect\n", encoding="utf-8")
    with pytest.raises(DataFormatError, match="expected header"):
        list(read_annotations(path))


def test_read_annotations_rejects_bad_aspect(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("Entry\tGO\tAspect\tEvidence\nP1\tGO:1\tX\tIDA\n", encoding="utf-8")
    with pytest.raises(DataFormatError, match="unknown aspect"):
        list(read_annotations(path))


def test_read_annotations_rejects_short_row(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("Entry\tGO\tAspect\tEvidence\nP1\tGO:1\n", encoding="utf-8")
    with pytest.raises(DataFormatError, match="expected 4 fields"):
        list(read_annotations(path))


def test_read_annotations_rejects_empty_file(tmp_path: Path) -> None:
    path = tmp_path / "empty.tsv"
    path.write_text("", encoding="utf-8")
    with pytest.raises(DataFormatError, match="empty annotation file"):
        list(read_annotations(path))


def test_load_labels_applies_curated_policy(fixtures_dir: Path) -> None:
    path = fixtures_dir / "human_annotations_mini.tsv"
    assert load_labels(path, aspect="F") == {"P1": frozenset({"GO:0000001"})}
    assert load_labels(path, aspect="P") == {"P2": frozenset({"GO:0000003"})}
    assert load_labels(path, aspect="C") == {}


def test_load_labels_rejects_bad_aspect(fixtures_dir: Path) -> None:
    with pytest.raises(DataFormatError, match="unknown aspect"):
        load_labels(fixtures_dir / "human_annotations_mini.tsv", aspect="X")
