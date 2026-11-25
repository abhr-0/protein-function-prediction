"""Repository layout and data paths.

The repository root is inferred from this file (``src/pfpred/paths.py``) and can
be overridden with ``PFP_ROOT`` for tests, Colab, or a non-editable install.
"""

from __future__ import annotations

import hashlib
import os
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

ENV_ROOT = "PFP_ROOT"
DEFAULT_ROOT = Path(__file__).resolve().parents[2]

GROUNDTRUTH_KINDS = ("NK", "LK", "PK", "PK_known")

HUMAN_FASTA_NAME = "human_proteins_2026_04.fasta.gz"
HUMAN_ANNOTATIONS_NAME = "human_go_annotations_2026_04.tsv.gz"


def repo_root(environ: Mapping[str, str] | None = None) -> Path:
    """Return the project root, honouring ``PFP_ROOT`` when it is set."""
    env = os.environ if environ is None else environ
    override = env.get(ENV_ROOT)
    return Path(override).expanduser().resolve() if override else DEFAULT_ROOT


@dataclass(frozen=True, slots=True)
class Layout:
    """Paths under a project root.

    Nothing here checks that files exist; adapters fail when they try to open a
    missing file.
    """

    root: Path

    @classmethod
    def from_env(cls, environ: Mapping[str, str] | None = None) -> Layout:
        return cls(repo_root(environ))

    @property
    def data(self) -> Path:
        return self.root / "data"

    @property
    def raw(self) -> Path:
        return self.data / "raw"

    @property
    def derived(self) -> Path:
        return self.data / "derived"

    @property
    def artifacts(self) -> Path:
        return self.root / "artifacts"

    @property
    def human_dir(self) -> Path:
        return self.raw / "uniprot"

    @property
    def human_fasta(self) -> Path:
        return self.human_dir / HUMAN_FASTA_NAME

    @property
    def human_annotations(self) -> Path:
        return self.human_dir / HUMAN_ANNOTATIONS_NAME

    @property
    def cafa6_dir(self) -> Path:
        return self.raw / "cafa6"

    @property
    def cafa6_train_terms(self) -> Path:
        return self.cafa6_dir / "train" / "train_terms.tsv"

    @property
    def cafa6_train_taxonomy(self) -> Path:
        return self.cafa6_dir / "train" / "train_taxonomy.tsv"

    @property
    def cafa6_train_sequences(self) -> Path:
        return self.cafa6_dir / "train" / "train_sequences.fasta"

    @property
    def cafa6_test_sequences(self) -> Path:
        return self.cafa6_dir / "eval" / "test_sequences.fasta"

    def cafa6_groundtruth(self, kind: str) -> Path:
        if kind not in GROUNDTRUTH_KINDS:
            raise ValueError(
                f"unknown ground-truth kind {kind!r}; expected one of {GROUNDTRUTH_KINDS}"
            )
        return self.cafa6_dir / "eval" / f"May2026_groundtruth_{kind}.tsv"


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    """Return the SHA-256 hex digest of ``path``."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()
