from __future__ import annotations

from pathlib import Path

import pytest

from pfpred.cafa6 import (
    groundtruth_proteins,
    groundtruth_union,
    read_groundtruth,
    read_sequences,
    read_taxonomy,
    read_terms,
)
from pfpred.records import Annotation, DataFormatError


def test_read_terms(fixtures_dir: Path) -> None:
    rows = list(read_terms(fixtures_dir / "cafa6_train_terms_mini.tsv"))
    assert rows[0] == Annotation(entry="P1", term="GO:0000001", aspect="F")
    assert len(rows) == 3
    assert all(row.evidence == "" for row in rows)


def test_read_terms_gzip(fixtures_dir: Path, gzipped) -> None:
    packed = gzipped(fixtures_dir / "cafa6_train_terms_mini.tsv")
    assert len(list(read_terms(packed))) == 3


def test_read_terms_rejects_wrong_header(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("EntryID\tterm\nP1\tGO:1\n", encoding="utf-8")
    with pytest.raises(DataFormatError, match="expected header"):
        list(read_terms(path))


def test_read_terms_rejects_bad_aspect(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("EntryID\tterm\taspect\nP1\tGO:1\tX\n", encoding="utf-8")
    with pytest.raises(DataFormatError, match="unknown aspect"):
        list(read_terms(path))


def test_read_taxonomy(fixtures_dir: Path) -> None:
    assert read_taxonomy(fixtures_dir / "cafa6_train_taxonomy_mini.tsv") == {
        "P1": 9606,
        "P2": 10090,
    }


def test_read_taxonomy_rejects_conflicting_taxids(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("P1\t9606\nP1\t10090\n", encoding="utf-8")
    with pytest.raises(DataFormatError, match="conflicting taxids"):
        read_taxonomy(path)


@pytest.mark.parametrize("taxid", ["abc", "-1", "0"])
def test_read_taxonomy_rejects_bad_taxid(tmp_path: Path, taxid: str) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text(f"P1\t{taxid}\n", encoding="utf-8")
    with pytest.raises(DataFormatError, match="taxid"):
        read_taxonomy(path)


def test_read_taxonomy_rejects_short_row(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("P1\n", encoding="utf-8")
    with pytest.raises(DataFormatError, match="expected 2 fields"):
        read_taxonomy(path)


def test_groundtruth_proteins(fixtures_dir: Path) -> None:
    assert groundtruth_proteins(fixtures_dir / "cafa6_groundtruth_PK_mini.tsv") == frozenset(
        {"P2", "P3"}
    )


def test_groundtruth_union_excludes_pk_known(fixtures_dir: Path) -> None:
    paths = [
        fixtures_dir / f"cafa6_groundtruth_{kind}_mini.tsv" for kind in ("NK", "LK", "PK")
    ]
    union = groundtruth_union(paths)
    assert union == frozenset({"P1", "P2", "P3"})
    with_known = [*paths, fixtures_dir / "cafa6_groundtruth_PK_known_mini.tsv"]
    assert groundtruth_union(with_known) == union


def test_read_groundtruth_kind(fixtures_dir: Path) -> None:
    rows = list(read_groundtruth(fixtures_dir / "cafa6_groundtruth_NK_mini.tsv"))
    assert rows == [Annotation(entry="P1", term="GO:0000001", aspect="F")]


def test_read_sequences_accepts_both_header_styles(fixtures_dir: Path) -> None:
    proteins = list(read_sequences(fixtures_dir / "cafa6_sequences_mini.fasta"))
    assert [(p.accession, p.sequence) for p in proteins] == [
        ("A0A0C5B5G6", "MRWQ"),
        ("A0A1B0GTW7", "MLLL"),
    ]
