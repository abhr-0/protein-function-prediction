from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from pfpred.splits import (
    SplitError,
    _is_windows_batch,
    _repair_and_check_leakage,
    _run_mmseqs,
    assign_cluster_splits,
    load_split,
    read_cluster_tsv,
    run_clustering,
    validate_split_payload,
)


def test_read_cluster_tsv_groups_members_and_singletons(tmp_path: Path) -> None:
    cluster_file = tmp_path / "clusters.tsv"
    cluster_file.write_text("P1\tP1\nP1\tP2\nP3\tP3\n", encoding="utf-8")
    assert read_cluster_tsv(cluster_file, {"P1", "P2", "P3"}) == {
        "P1": ("P1", "P2"), "P3": ("P3",),
    }


def test_read_cluster_tsv_rejects_duplicate_members(tmp_path: Path) -> None:
    path = tmp_path / "clusters.tsv"
    path.write_text("P1\tP1\nP1\tP2\nP3\tP2\n", encoding="utf-8")
    with pytest.raises(SplitError, match="duplicate cluster member"):
        read_cluster_tsv(path, {"P1", "P2", "P3"})


def test_read_cluster_tsv_rejects_missing_and_unknown_proteins(tmp_path: Path) -> None:
    path = tmp_path / "clusters.tsv"
    path.write_text("P1\tP1\nP1\tP9\n", encoding="utf-8")
    with pytest.raises(SplitError, match="unknown member"):
        read_cluster_tsv(path, {"P1", "P2"})
    path.write_text("P1\tP1\n", encoding="utf-8")
    with pytest.raises(SplitError, match="missing from cluster output"):
        read_cluster_tsv(path, {"P1", "P2"})


def test_read_cluster_tsv_rejects_wrong_columns_and_empty_file(tmp_path: Path) -> None:
    path = tmp_path / "clusters.tsv"
    path.write_text("P1\tP1\textra\n", encoding="utf-8")
    with pytest.raises(SplitError, match="expected representative"):
        read_cluster_tsv(path, {"P1"})
    path.write_text("", encoding="utf-8")
    with pytest.raises(SplitError, match="empty cluster file"):
        read_cluster_tsv(path, {"P1"})


def test_read_cluster_tsv_rejects_representative_without_self_membership(
    tmp_path: Path,
) -> None:
    path = tmp_path / "clusters.tsv"
    path.write_text("P1\tP2\nP3\tP3\n", encoding="utf-8")
    with pytest.raises(SplitError, match="absent from its own cluster"):
        read_cluster_tsv(path, {"P1", "P2", "P3"})


def test_assign_cluster_splits_is_order_independent_and_keeps_clusters_whole() -> None:
    clusters = {
        "P1": ("P1", "P2", "P3"), "P4": ("P4", "P5"),
        "P6": ("P6",), "P7": ("P7",), "P8": ("P8",), "P9": ("P9",),
    }
    first = assign_cluster_splits(clusters, seed=17)
    second = assign_cluster_splits(dict(reversed(list(clusters.items()))), seed=17)
    assert first == second
    assert set(first) == set(clusters)


def test_assign_cluster_splits_rejects_overlap_empty_and_too_few_clusters() -> None:
    with pytest.raises(SplitError, match="empty cluster collection"):
        assign_cluster_splits({})
    with pytest.raises(SplitError, match="more than one cluster"):
        assign_cluster_splits({"P1": ("P1", "P2"), "P2": ("P2", "P3")})
    with pytest.raises(SplitError, match="too small"):
        assign_cluster_splits({"P1": ("P1",), "P2": ("P2",)})


def _valid_payload() -> dict[str, object]:
    assignments = {
        "P1": {"cluster": "P1", "split": "train"},
        "P2": {"cluster": "P1", "split": "train"},
        "P3": {"cluster": "P3", "split": "validation"},
        "P4": {"cluster": "P4", "split": "test"},
    }
    return {
        "schema_version": 2,
        "input_sha256": "digest",
        "split": {
            "protein_counts": {"train": 2, "validation": 1, "test": 1},
            "cluster_counts": {"train": 1, "validation": 1, "test": 1},
            "total_proteins": 4,
            "total_clusters": 3,
        },
        "leak_check": {"hits": 0},
        "assignments": assignments,
    }


def test_validate_split_payload_accepts_consistent_payload() -> None:
    validate_split_payload(_valid_payload())


def test_validate_split_payload_rejects_split_cluster_and_count_mismatch() -> None:
    payload = _valid_payload()
    payload["assignments"]["P2"]["split"] = "test"
    with pytest.raises(SplitError, match="divided across splits"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["split"]["total_proteins"] = 99
    with pytest.raises(SplitError, match="total_proteins"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["hits"] = 2
    with pytest.raises(SplitError, match="zero-hit"):
        validate_split_payload(payload)


def test_load_split_checks_fasta_digest(tmp_path: Path) -> None:
    payload = _valid_payload()
    path = tmp_path / "split.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">P1\nAAA\n", encoding="utf-8")
    with pytest.raises(SplitError, match="input_sha256"):
        load_split(path, fasta_path=fasta)


def test_run_clustering_writes_provenance_and_zero_hit_check(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    fasta = tmp_path / "mini.fasta"
    fasta.write_text(">P1\nACD\n>P2\nACE\n>P3\nFGH\n", encoding="utf-8")
    output = tmp_path / "derived" / "human_split.json"

    def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        if command[1] == "version":
            return SimpleNamespace(stdout="MMseqs2 test-version", stderr="")
        if command[1] == "easy-cluster":
            assert "--cluster-mode" in command
            assert "--cluster-reassign" in command
            prefix = Path(command[3])
            Path(f"{prefix}_cluster.tsv").write_text(
                "P1\tP1\nP2\tP2\nP3\tP3\n", encoding="utf-8",
            )
        elif command[1] == "easy-search":
            Path(command[4]).write_text("", encoding="utf-8")
        else:
            raise AssertionError(command)
        return SimpleNamespace(stdout="", stderr="")

    monkeypatch.setattr("pfpred.splits.subprocess.run", fake_run)
    payload = run_clustering(fasta, output, seed=11)
    assert json.loads(output.read_text(encoding="utf-8")) == payload
    assert payload["mmseqs"]["version"] == "MMseqs2 test-version"
    assert payload["mmseqs"]["parameters"]["cluster_mode"] == 0
    assert payload["mmseqs"]["parameters"]["cluster_reassign"] is True
    assert payload["leak_check"]["hits"] == 0
    assert payload["split"]["total_proteins"] == 3
    assert payload["split"]["total_clusters"] == 3
    assert load_split(output, fasta_path=fasta) == payload


def test_windows_batch_detection_and_invocation(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr("pfpred.splits._is_windows_batch", lambda _: True)
    monkeypatch.setattr(
        "pfpred.splits.subprocess.run",
        lambda command, **kwargs: calls.append(command) or SimpleNamespace(stdout="", stderr=""),
    )
    _run_mmseqs(r"C:\Program Files (x86)\mmseqs\mmseqs.bat", ["version"], check=True)
    assert calls == [[
        "cmd.exe", "/d", "/c", "call",
        r"C:\Program Files (x86)\mmseqs\mmseqs.bat", "version",
    ]]



def test_leakage_repair_moves_entire_offending_cluster_to_train(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    clusters = {
        "P1": ("P1",),
        "P2": ("P2", "P2b"),
        "P2c": ("P2c",),
        "P3": ("P3",),
    }
    assignments = {"P1": "train", "P2": "validation", "P2c": "validation", "P3": "test"}
    sequences = {name: "ACDE" for members in clusters.values() for name in members}
    calls = 0

    def fake_search(
        executable: str, query: Path, target: Path, work_dir: Path, label: str,
    ) -> set[str]:
        nonlocal calls
        calls += 1
        return {"P2"} if calls == 1 else set()

    monkeypatch.setattr("pfpred.splits._search_hits", fake_search)
    result = _repair_and_check_leakage(
        "mmseqs", clusters, assignments, sequences, tmp_path,
    )

    assert assignments == {"P1": "train", "P2": "train", "P2c": "validation", "P3": "test"}
    assert result["hits"] == 0
    assert result["moved_clusters_to_train"] == 1
    assert result["moved_proteins_to_train"] == 2


def test_validate_split_payload_rejects_empty_partition() -> None:
    payload = _valid_payload()
    payload["assignments"].pop("P4")
    payload["split"]["protein_counts"]["test"] = 0
    payload["split"]["cluster_counts"]["test"] = 0
    payload["split"]["total_proteins"] = 3
    payload["split"]["total_clusters"] = 2
    with pytest.raises(SplitError, match="empty split"):
        validate_split_payload(payload)
