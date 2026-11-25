"""Shared test helpers."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures_dir() -> Path:
    return FIXTURES


@pytest.fixture
def gzipped(tmp_path: Path):
    """Return a callable that gzips a fixture into ``tmp_path``."""

    def _gzip(src: Path) -> Path:
        dst = tmp_path / f"{src.name}.gz"
        with open(src, "rb") as raw, gzip.open(dst, "wb") as packed:
            shutil.copyfileobj(raw, packed)
        return dst

    return _gzip
