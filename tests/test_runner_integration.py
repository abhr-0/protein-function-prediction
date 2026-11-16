"""End-to-end wiring test for the baseline runner on synthetic data.

Exercises ``--fasta-path`` + ``--annotations-path`` against a tiny generated
UniProt-style FASTA and flat file, so the test never touches real data and
never writes into ``data/artifacts``.
"""

import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_baseline_subset.py"

AA = "ACDEFGHIKLMNPQRSTVWY"


@pytest.fixture(scope="module")
def runner():
    """Load scripts/run_baseline_subset.py as a module."""
    spec = importlib.util.spec_from_file_location("run_baseline_subset_it", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _write_fasta(path: Path, n_proteins: int = 24) -> None:
    lines = []
    for index in range(n_proteins):
        accession = f"Q{index:05d}"
        sequence = (AA * (1 + index % 4))[: 20 + index]
        lines.append(
            f">sp|{accession}|PROT{index}_HUMAN Protein {index} "
            f"OS=Homo sapiens OX=9606 PE=1 SV=1\n{sequence}\n"
        )
    path.write_text("".join(lines))


def _write_dat(path: Path, n_proteins: int = 24) -> None:
    """Three shared terms plus one unique term per protein."""
    entries = []
    for index in range(n_proteins):
        accession = f"Q{index:05d}"
        go_lines = [
            "DR   GO; GO:0000001; F:function one; IDA:UniProtKB.",
            "DR   GO; GO:0000002; F:function two; IDA:UniProtKB.",
        ]
        if index % 3 == 0:
            go_lines.append("DR   GO; GO:0000003; F:function three; IDA:UniProtKB.")
        go_lines.append(f"DR   GO; GO:9{index:06d}; F:unique {index}; IDA:UniProtKB.")
        entries.append(
            f"ID   PROT{index}_HUMAN             Reviewed;         100 AA.\n"
            f"AC   {accession};\n" + "\n".join(go_lines) + "\n//\n"
        )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("".join(entries))


@pytest.fixture()
def synthetic_inputs(tmp_path: Path) -> tuple[Path, Path]:
    fasta = tmp_path / "synthetic.fasta"
    dat = tmp_path / "synthetic.dat.gz"
    _write_fasta(fasta)
    _write_dat(dat)
    return fasta, dat


def test_runner_runs_end_to_end_on_synthetic_inputs(
    runner, synthetic_inputs, tmp_path, monkeypatch
):
    """The runner must wire --fasta-path/--annotations-path through to the model."""
    fasta, dat = synthetic_inputs
    artifacts = tmp_path / "artifacts"
    artifacts.mkdir()

    # Keep every artifact, cache, and model inside tmp_path. Patch the names
    # the script itself bound, otherwise the annotation cache escapes to the
    # real artifacts directory.
    cache_path = tmp_path / "cache" / "annotations.tsv.gz"
    monkeypatch.setattr(runner, "get_artifacts_dir", lambda: artifacts)
    monkeypatch.setattr(runner, "get_annotations_cache_path", lambda **_: cache_path)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_baseline_subset.py",
            "--fasta-path", str(fasta),
            "--annotations-path", str(dat),
            "--n-proteins", "100",
            "--min-term-count", "2",
            "--top-terms", "0",
            "--test-size", "0.25",
            "--val-size", "0.1",
            "--seed", "7",
            "--aspect", "F",
        ],
    )

    runner.main()

    metrics_files = list(artifacts.glob("baseline_metrics_*.json"))
    assert len(metrics_files) == 1
    metrics = json.loads(metrics_files[0].read_text())

    assert metrics["config"]["fasta_file"] == str(fasta)
    assert metrics["config"]["label_file"] == str(dat)
    assert metrics["config"]["seed"] == 7
    assert metrics["counts"]["matched_proteins"] == 24
    assert metrics["counts"]["target_terms"] >= 2
    assert metrics["counts"]["test_proteins"] > 0
    assert metrics["threshold"]["at_search_floor"] in (True, False)
    # Both artifacts share the run tag, so they stay associated.
    assert list(artifacts.glob("baseline_model_*.joblib"))

    # The unique per-protein terms must not have inflated the vocabulary.
    assert metrics["counts"]["target_terms"] < 10


def test_runner_rejects_missing_fasta(runner, tmp_path, monkeypatch):
    """A missing sequence file must fail with an actionable message."""
    dat = tmp_path / "synthetic.dat.gz"
    _write_dat(dat)
    monkeypatch.setattr(
        runner, "get_annotations_cache_path", lambda **_: tmp_path / "cache.tsv.gz"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_baseline_subset.py",
            "--fasta-path", str(tmp_path / "absent.fasta"),
            "--annotations-path", str(dat),
            "--n-proteins", "10",
        ],
    )

    with pytest.raises(FileNotFoundError, match="Sequence file not found"):
        runner.main()
