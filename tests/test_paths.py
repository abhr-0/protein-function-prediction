from __future__ import annotations

from pathlib import Path

import pytest

from pfpred.paths import DEFAULT_ROOT, Layout, repo_root, sha256_file


def test_default_root_contains_pyproject() -> None:
    assert (DEFAULT_ROOT / "pyproject.toml").is_file()


def test_repo_root_honours_override(tmp_path: Path) -> None:
    assert repo_root({"PFP_ROOT": str(tmp_path)}) == tmp_path.resolve()


def test_layout_composes_paths(tmp_path: Path) -> None:
    layout = Layout(tmp_path)
    assert layout.raw == tmp_path / "data" / "raw"
    assert layout.derived == tmp_path / "data" / "derived"
    assert layout.artifacts == tmp_path / "artifacts"
    assert layout.human_fasta.name == "human_proteins_2026_04.fasta.gz"
    assert layout.human_annotations.name == "human_go_annotations_2026_04.tsv.gz"
    assert layout.cafa6_train_terms == tmp_path / "data/raw/cafa6/train/train_terms.tsv"
    assert layout.cafa6_test_sequences == tmp_path / "data/raw/cafa6/eval/test_sequences.fasta"


@pytest.mark.parametrize("kind", ["NK", "LK", "PK", "PK_known"])
def test_groundtruth_paths(tmp_path: Path, kind: str) -> None:
    path = Layout(tmp_path).cafa6_groundtruth(kind)
    assert path.name == f"May2026_groundtruth_{kind}.tsv"


def test_unknown_groundtruth_kind_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unknown ground-truth kind"):
        Layout(tmp_path).cafa6_groundtruth("XX")


def test_sha256_file(tmp_path: Path) -> None:
    target = tmp_path / "abc.txt"
    target.write_bytes(b"abc")
    assert sha256_file(target) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
