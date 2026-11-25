"""Shared test helpers."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


def _gzip_copy(src: Path, dst: Path) -> None:
    with open(src, "rb") as raw, gzip.open(dst, "wb") as packed:
        shutil.copyfileobj(raw, packed)


@pytest.fixture
def fixtures_dir() -> Path:
    return FIXTURES


@pytest.fixture
def gzipped(tmp_path: Path):
    """Return a callable that gzips a fixture into ``tmp_path``."""

    def _gzip(src: Path) -> Path:
        dst = tmp_path / f"{src.name}.gz"
        _gzip_copy(src, dst)
        return dst

    return _gzip


@pytest.fixture
def mini_layout(tmp_path: Path) -> Path:
    """Build a miniature project tree with the expected file names."""
    human = tmp_path / "data" / "raw" / "uniprot"
    human.mkdir(parents=True)
    _gzip_copy(FIXTURES / "human_proteins_mini.fasta", human / "human_proteins_2026_04.fasta.gz")
    _gzip_copy(
        FIXTURES / "human_annotations_mini.tsv",
        human / "human_go_annotations_2026_04.tsv.gz",
    )

    train = tmp_path / "data" / "raw" / "cafa6" / "train"
    train.mkdir(parents=True)
    shutil.copy(FIXTURES / "cafa6_train_terms_mini.tsv", train / "train_terms.tsv")
    shutil.copy(FIXTURES / "cafa6_train_taxonomy_mini.tsv", train / "train_taxonomy.tsv")

    evaluation = tmp_path / "data" / "raw" / "cafa6" / "eval"
    evaluation.mkdir(parents=True)
    for kind in ("NK", "LK", "PK", "PK_known"):
        shutil.copy(
            FIXTURES / f"cafa6_groundtruth_{kind}_mini.tsv",
            evaluation / f"May2026_groundtruth_{kind}.tsv",
        )
    shutil.copy(FIXTURES / "cafa6_sequences_mini.fasta", evaluation / "test_sequences.fasta")
    return tmp_path
