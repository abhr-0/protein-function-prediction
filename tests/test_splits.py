from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from pfpred.splits import (
    _run_mmseqs,
    SplitError,
    assign_cluster_splits,
    read_cluster_tsv,
    run_clustering,
)


def test_read_cluster_tsv_groups_members_and_singletons(tmp_path: Path) -> None:
    cluster_file = tmp_path / "clusters.tsv"
    cluster_file.write_text(
        "P1\tP1\nP1\tP2\nP3\tP3\n",
        encoding="utf-8",
    )

    assert read_cluster_tsv(cluster_file, {"P1", "P2", "P3"}) == {
        "P1": ("P1", "P2"),
        "P3": ("P3",),
    }


def test_read_cluster_tsv_rejects_duplicate_members(tmp_path: Path) -> None:
    cluster_file = tmp_path / "clusters.tsv"
    cluster_file.write_text("P1\tP1\nP1\tP2\nP3\tP2\n", encoding="utf-8")

    with pytest.raises(SplitError, match="duplicate cluster member"):
        read_cluster_tsv(cluster_file, {"P1", "P2", "P3"})


def test_read_cluster_tsv_rejects_missing_input_proteins(tmp_path: Path) -> None:
    cluster_file = tmp_path / "clusters.tsv"
    cluster_file.write_text("P1\tP1\n", encoding="utf-8")

    with pytest.raises(SplitError, match="missing from cluster output"):
        read_cluster_tsv(cluster_file, {"P1", "P2"})


def test_read_cluster_tsv_rejects_unknown_proteins(tmp_path: Path) -> None:
    cluster_file = tmp_path / "clusters.tsv"
    cluster_file.write_text("P1\tP1\nP1\tP9\n", encoding="utf-8")

    with pytest.raises(SplitError, match="unknown member"):
        read_cluster_tsv(cluster_file, {"P1", "P2"})


def test_read_cluster_tsv_rejects_representative_without_self_membership(
    tmp_path: Path,
) -> None:
    cluster_file = tmp_path / "clusters.tsv"
    cluster_file.write_text("P1\tP2\nP3\tP3\n", encoding="utf-8")

    with pytest.raises(SplitError, match="absent from its own cluster"):
        read_cluster_tsv(cluster_file, {"P1", "P2", "P3"})


def test_assign_cluster_splits_is_deterministic_and_keeps_clusters_whole() -> None:
    clusters = {
        "P1": ("P1", "P2", "P3"),
        "P4": ("P4", "P5"),
        "P6": ("P6",),
        "P7": ("P7",),
        "P8": ("P8",),
        "P9": ("P9",),
    }

    first = assign_cluster_splits(clusters, seed=17)
    second = assign_cluster_splits(clusters, seed=17)

    assert first == second
    assert set(first) == set(clusters)
    protein_splits = {
        member: first[representative]
        for representative, members in clusters.items()
        for member in members
    }
    for members in clusters.values():
        assert len({protein_splits[member] for member in members}) == 1


def test_assign_cluster_splits_rejects_overlapping_clusters() -> None:
    clusters = {"P1": ("P1", "P2"), "P2": ("P2", "P3")}

    with pytest.raises(SplitError, match="more than one cluster"):
        assign_cluster_splits(clusters)


def test_assign_cluster_splits_rejects_empty_input() -> None:
    with pytest.raises(SplitError, match="empty cluster collection"):
        assign_cluster_splits({})


def test_run_clustering_writes_provenance_and_assignments(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fasta = tmp_path / "mini.fasta"
    fasta.write_text(">P1\nACD\n>P2\nACE\n>P3\nFGH\n", encoding="utf-8")
    output = tmp_path / "derived" / "human_split.json"

    def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        if command[1] == "version":
            return SimpleNamespace(stdout="MMseqs2 test-version", stderr="")
        assert command[1] == "easy-cluster"
        assert command[6:] == ["0.30", "-c", "0.80", "--cov-mode", "0"]
        prefix = Path(command[3])
        Path(f"{prefix}_cluster.tsv").write_text(
            "P1\tP1\nP1\tP2\nP3\tP3\n",
            encoding="utf-8",
        )
        return SimpleNamespace(stdout="", stderr="")

    monkeypatch.setattr("pfpred.splits.subprocess.run", fake_run)

    payload = run_clustering(fasta, output, seed=11)

    assert output.exists()
    assert json.loads(output.read_text(encoding="utf-8")) == payload
    assert payload["mmseqs"]["version"] == "MMseqs2 test-version"
    assert payload["mmseqs"]["parameters"] == {
        "min_seq_id": 0.30,
        "coverage": 0.80,
        "coverage_mode": 0,
    }
    assert payload["split"]["total_proteins"] == 3
    assert payload["split"]["total_clusters"] == 2



def test_run_mmseqs_uses_cmd_for_windows_batch_files(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[list[str]] = []

    def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        calls.append(command)
        return SimpleNamespace(stdout="MMseqs2 test-version", stderr="")

    monkeypatch.setattr("pfpred.splits.os.name", "nt")
    monkeypatch.setattr("pfpred.splits.subprocess.run", fake_run)

    _run_mmseqs(
        r"C:\Program Files (x86)\mmseqs\mmseqs.bat",
        ["version"],
        check=True,
        capture_output=True,
        text=True,
    )

    assert calls == [
        [
            "cmd.exe",
            "/d",
            "/c",
            "call",
            r"C:\Program Files (x86)\mmseqs\mmseqs.bat",
            "version",
        ]
    ]
