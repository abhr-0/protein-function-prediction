from __future__ import annotations

import json
from pathlib import Path

import pytest

from pfpred.cli import EXPECTED, collect_counts, compare, main

MINI_COUNTS = {
    "human.proteins": 3,
    "human.residues": 12,
    "human.annotations": 3,
    "human.terms": 3,
    "human.proteins_with_labels": 2,
    "human.excluded_rows": 1,
    "human.curated.annotations": 2,
    "human.curated.proteins": 2,
    "human.curated.terms": 2,
    "human.curated.F": 1,
    "human.curated.P": 1,
    "human.curated.C": 0,
    "cafa6.train.proteins": 2,
    "cafa6.train.annotations": 3,
    "cafa6.train.F": 1,
    "cafa6.train.P": 1,
    "cafa6.train.C": 1,
    "cafa6.species": 2,
    "gt.NK": 1,
    "gt.LK": 1,
    "gt.PK": 2,
    "gt.union": 3,
    "gt.LK_int_PK": 1,
    "gt.in_train": 2,
    "cafa6.test_sequences": 2,
}


def test_collect_counts_on_fixture_tree(mini_layout: Path) -> None:
    from pfpred.paths import Layout

    assert collect_counts(Layout(mini_layout)) == MINI_COUNTS


def test_compare_accepts_identical() -> None:
    assert compare(MINI_COUNTS, MINI_COUNTS) == []


def test_compare_reports_drift_missing_and_unexpected() -> None:
    actual = {**MINI_COUNTS, "human.proteins": 4}
    del actual["human.terms"]
    actual["extra.key"] = 1
    problems = compare(actual, MINI_COUNTS)
    assert "human.proteins: expected 3, got 4" in problems
    assert "human.terms: missing" in problems
    assert "extra.key: unexpected" in problems


def test_main_no_check_prints_table(mini_layout: Path, capsys) -> None:
    assert main(["verify", "--root", str(mini_layout), "--no-check"]) == 0
    out = capsys.readouterr().out
    assert "human.proteins" in out
    assert "3" in out


def test_main_json(mini_layout: Path, capsys) -> None:
    assert main(["verify", "--root", str(mini_layout), "--json", "--no-check"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload == MINI_COUNTS


def test_main_flags_drift(mini_layout: Path, capsys) -> None:
    assert main(["verify", "--root", str(mini_layout)]) == 1
    err = capsys.readouterr().err
    assert "error: human.proteins: expected 20429, got 3" in err


def test_main_requires_a_command() -> None:
    with pytest.raises(SystemExit):
        main([])


def test_expected_has_no_nested_mismatch() -> None:
    assert set(EXPECTED) >= {"human.proteins", "cafa6.train.proteins", "gt.union"}


def test_cluster_command_invokes_clustering(
    mini_layout: Path, monkeypatch: pytest.MonkeyPatch, capsys
) -> None:
    from pfpred import cli
    from pfpred.paths import Layout

    captured: dict[str, object] = {}

    def fake_run_clustering(
        fasta_path: Path,
        output_path: Path,
        *,
        mmseqs_executable: str,
        seed: int,
    ) -> dict[str, object]:
        captured.update(
            fasta_path=fasta_path,
            output_path=output_path,
            mmseqs_executable=mmseqs_executable,
            seed=seed,
        )
        return {"split": {"total_proteins": 3}}

    monkeypatch.setattr(cli, "run_clustering", fake_run_clustering)

    assert main(["cluster", "--root", str(mini_layout), "--seed", "9"]) == 0
    layout = Layout(mini_layout)
    assert captured == {
        "fasta_path": layout.human_fasta,
        "output_path": layout.derived / "human_split.json",
        "mmseqs_executable": "mmseqs",
        "seed": 9,
    }
    assert "total_proteins" in capsys.readouterr().out



def test_cluster_command_reports_malformed_fasta_without_traceback(
    mini_layout: Path, capsys,
) -> None:
    from pfpred.paths import Layout

    layout = Layout(mini_layout)
    layout.human_fasta.parent.mkdir(parents=True, exist_ok=True)
    layout.human_fasta.write_text("not a FASTA record\n", encoding="utf-8")

    assert main(["cluster", "--root", str(mini_layout)]) == 1
    captured = capsys.readouterr()
    assert "clustering failed" in captured.err
    assert "Traceback" not in captured.err
