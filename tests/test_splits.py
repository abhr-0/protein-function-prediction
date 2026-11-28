from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from pfpred.splits import (
    SplitError,
    _canonical_executable,
    _repair_and_check_leakage,
    _run_mmseqs,
    _search_hits,
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


def _mmseqs_result_writer(outputs: list[str]):
    """Fake MMseqs2 easy-search that writes the next canned result file."""
    pending = iter(outputs)

    def fake_run(executable: str, arguments: list[str], **kwargs: object) -> SimpleNamespace:
        Path(arguments[3]).write_text(next(pending), encoding="utf-8")
        return SimpleNamespace(stdout="", stderr="")

    return fake_run


def test_search_hits_returns_query_identity_pairs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "pfpred.splits._run_mmseqs",
        _mmseqs_result_writer(["P1\tP2\t0.875\nP1\tP3\t0.5\n"]),
    )
    assert _search_hits("mmseqs", tmp_path / "q.fasta", tmp_path / "t.fasta", tmp_path, "ok") == [
        ("P1", 0.875), ("P1", 0.5),
    ]


def test_search_hits_rejects_malformed_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "pfpred.splits._run_mmseqs",
        _mmseqs_result_writer([
            "P1\tP2\n", "P1\tP2\tnot-a-number\n", "P1\t\t0.5\n",
            "P1\tP2\tnan\n", "P1\tP2\tinf\n",
        ]),
    )
    for label, match in (
        ("columns", "malformed"),
        ("identity", "non-numeric identity"),
        ("empty", "malformed"),
        ("nan", "non-finite identity"),
        ("inf", "non-finite identity"),
    ):
        with pytest.raises(SplitError, match=match):
            _search_hits("mmseqs", tmp_path / "q.fasta", tmp_path / "t.fasta", tmp_path, label)


def test_search_hits_rejects_duplicate_label(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr("pfpred.splits._run_mmseqs", _mmseqs_result_writer(["", ""]))
    _search_hits("mmseqs", tmp_path / "q.fasta", tmp_path / "t.fasta", tmp_path, "dup")
    with pytest.raises(SplitError, match="duplicate"):
        _search_hits("mmseqs", tmp_path / "q.fasta", tmp_path / "t.fasta", tmp_path, "dup")


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


def test_assign_cluster_splits_breaks_zero_fill_ties_toward_train() -> None:
    clusters = {
        "P1": ("P1", "P2", "P3", "P4"),
        "P5": ("P5",),
        "P6": ("P6",),
        "P7": ("P7",),
    }
    assert assign_cluster_splits(clusters, seed=3)["P1"] == "train"


def _valid_payload() -> dict[str, object]:
    assignments = {
        "P1": {"cluster": "P1", "split": "train"},
        "P2": {"cluster": "P1", "split": "train"},
        "P3": {"cluster": "P3", "split": "validation"},
        "P4": {"cluster": "P4", "split": "test"},
    }
    gate = {
        "min_seq_id": 0.30,
        "coverage": 0.80,
        "coverage_mode": 0,
        "max_seqs": 1000,
    }
    return {
        "schema_version": 3,
        "input_sha256": "a" * 64,
        "platform": {"system": "linux", "machine": "x86_64", "python": "3.12.0"},
        "mmseqs": {
            "executable": "/usr/bin/mmseqs",
            "version": "18-8cc5c",
            "workflow": "easy-cluster",
            "parameters": {
                "min_seq_id": 0.30,
                "coverage": 0.80,
                "coverage_mode": 0,
                "cluster_mode": 0,
                "cluster_reassign": True,
                "search_sensitivity": 5.0,
                "search_max_seqs": 1000,
            },
        },
        "split": {
            "protein_counts": {"train": 2, "validation": 1, "test": 1},
            "cluster_counts": {"train": 1, "validation": 1, "test": 1},
            "total_proteins": 4,
            "total_clusters": 3,
        },
        "leak_check": {
            "hits": 0,
            **gate,
            "sensitivity": 5.0,
            "rounds": 1,
            "hits_observed_before_repair": 0,
            "moved_clusters_to_train": 0,
            "moved_proteins_to_train": 0,
            "final_protein_counts": {"train": 2, "validation": 1, "test": 1},
            "holdout_cross_check": {
                **gate,
                "sensitivity": 5.0,
                "repaired": False,
                "validation_vs_test": {"queries_with_hits": 0, "max_identity": 0.0},
                "test_vs_validation": {"queries_with_hits": 0, "max_identity": 0.0},
            },
            "strict_probe": {
                **gate,
                "sensitivity": 7.5,
                "repaired": False,
                "validation_vs_train": {"queries_with_hits": 0, "max_identity": 0.0},
                "test_vs_train": {"queries_with_hits": 0, "max_identity": 0.0},
            },
        },
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


def test_validate_split_payload_rejects_bad_schema_and_leak_sections() -> None:
    payload = _valid_payload()
    payload["schema_version"] = 2
    with pytest.raises(SplitError, match="schema_version"):
        validate_split_payload(payload)
    payload = _valid_payload()
    del payload["leak_check"]["holdout_cross_check"]
    with pytest.raises(SplitError, match="holdout_cross_check"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["strict_probe"]["test_vs_train"].pop("max_identity")
    with pytest.raises(SplitError, match="max_identity"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["strict_probe"]["validation_vs_train"]["queries_with_hits"] = "0"
    with pytest.raises(SplitError, match="queries_with_hits"):
        validate_split_payload(payload)


def test_validate_split_payload_rejects_unpinned_probe_parameters() -> None:
    payload = _valid_payload()
    payload["leak_check"]["hits"] = False
    with pytest.raises(SplitError, match="zero-hit"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["sensitivity"] = 7.5
    with pytest.raises(SplitError, match="leak_check.sensitivity"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["holdout_cross_check"]["min_seq_id"] = 0.5
    with pytest.raises(SplitError, match="min_seq_id"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["holdout_cross_check"]["repaired"] = True
    with pytest.raises(SplitError, match="repaired"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["strict_probe"]["sensitivity"] = 5.0
    with pytest.raises(SplitError, match="strict_probe.sensitivity"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["strict_probe"]["max_seqs"] = 2000
    with pytest.raises(SplitError, match="max_seqs"):
        validate_split_payload(payload)


def test_validate_split_payload_pins_coverage_mode_everywhere() -> None:
    for path in (
        ("leak_check",),
        ("leak_check", "holdout_cross_check"),
        ("leak_check", "strict_probe"),
    ):
        payload = _valid_payload()
        section = payload
        for key in path:
            section = section[key]
        section["coverage_mode"] = 1
        with pytest.raises(SplitError, match="coverage_mode"):
            validate_split_payload(payload)


def test_validate_split_payload_requires_cluster_definition() -> None:
    payload = _valid_payload()
    payload["mmseqs"]["version"] = ""
    with pytest.raises(SplitError, match="mmseqs.version"):
        validate_split_payload(payload)
    payload = _valid_payload()
    del payload["mmseqs"]["version"]
    with pytest.raises(SplitError, match="mmseqs.version"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["mmseqs"]["workflow"] = "easy-search"
    with pytest.raises(SplitError, match="workflow"):
        validate_split_payload(payload)
    payload = _valid_payload()
    del payload["mmseqs"]["parameters"]
    with pytest.raises(SplitError, match="parameters"):
        validate_split_payload(payload)
    for key, value in (
        ("cluster_mode", 1),
        ("cluster_reassign", False),
        ("search_sensitivity", 7.5),
        ("coverage_mode", 1),
    ):
        payload = _valid_payload()
        payload["mmseqs"]["parameters"][key] = value
        with pytest.raises(SplitError, match=key):
            validate_split_payload(payload)


def test_validate_split_payload_rejects_bad_leak_bookkeeping() -> None:
    payload = _valid_payload()
    payload["leak_check"]["rounds"] = 0
    with pytest.raises(SplitError, match="rounds"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["rounds"] = True
    with pytest.raises(SplitError, match="rounds must be an integer"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["moved_clusters_to_train"] = -1
    with pytest.raises(SplitError, match="moved_clusters_to_train"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["hits_observed_before_repair"] = 1.5
    with pytest.raises(SplitError, match="hits_observed_before_repair must be an integer"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["final_protein_counts"]["test"] = 2
    with pytest.raises(SplitError, match="final_protein_counts"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["final_protein_counts"].pop("test")
    with pytest.raises(SplitError, match="final_protein_counts"):
        validate_split_payload(payload)


def test_validate_summary_rejects_more_hits_than_split_proteins() -> None:
    for name, direction, query_split in (
        ("holdout_cross_check", "validation_vs_test", "validation"),
        ("holdout_cross_check", "test_vs_validation", "test"),
        ("strict_probe", "validation_vs_train", "validation"),
        ("strict_probe", "test_vs_train", "test"),
    ):
        payload = _valid_payload()
        summary = payload["leak_check"][name][direction]
        summary["queries_with_hits"] = 2
        summary["max_identity"] = 0.5
        with pytest.raises(SplitError, match=f"exceeds the {query_split} protein count"):
            validate_split_payload(payload)


def test_canonical_executable_resolves_symlinks_and_falls_back(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "store" / "mmseqs"
    target.parent.mkdir()
    target.write_text("", encoding="utf-8")
    link = tmp_path / "profile" / "bin" / "mmseqs"
    link.parent.mkdir(parents=True)
    try:
        link.symlink_to(target)
    except (NotImplementedError, OSError) as exc:
        pytest.skip(f"symlink creation is unavailable: {exc}")
    monkeypatch.setattr("pfpred.splits.shutil.which", lambda _name: str(link))
    assert _canonical_executable("mmseqs") == str(target.resolve())
    monkeypatch.setattr("pfpred.splits.shutil.which", lambda _name: None)
    assert _canonical_executable("mmseqs") == "mmseqs"


def test_validate_split_payload_requires_digest_platform_and_executable() -> None:
    payload = _valid_payload()
    payload["input_sha256"] = "digest"
    with pytest.raises(SplitError, match="input_sha256"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["input_sha256"] = "A" * 64
    with pytest.raises(SplitError, match="input_sha256"):
        validate_split_payload(payload)
    payload = _valid_payload()
    del payload["input_sha256"]
    with pytest.raises(SplitError, match="input_sha256"):
        validate_split_payload(payload)
    payload = _valid_payload()
    del payload["platform"]
    with pytest.raises(SplitError, match="platform"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["platform"]["machine"] = ""
    with pytest.raises(SplitError, match="platform"):
        validate_split_payload(payload)
    payload = _valid_payload()
    del payload["mmseqs"]["executable"]
    with pytest.raises(SplitError, match="executable"):
        validate_split_payload(payload)


def test_validate_split_payload_rejects_bad_identity_summaries() -> None:
    payload = _valid_payload()
    payload["leak_check"]["strict_probe"]["test_vs_train"]["queries_with_hits"] = -1
    with pytest.raises(SplitError, match="non-negative"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["holdout_cross_check"]["validation_vs_test"]["queries_with_hits"] = True
    with pytest.raises(SplitError, match="integer queries_with_hits"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["holdout_cross_check"]["validation_vs_test"]["max_identity"] = 1.5
    with pytest.raises(SplitError, match=r"within \[0\.0, 1\.0\]"):
        validate_split_payload(payload)
    payload = _valid_payload()
    summary = payload["leak_check"]["strict_probe"]["test_vs_train"]
    summary["queries_with_hits"] = 1
    summary["max_identity"] = float("nan")
    with pytest.raises(SplitError, match="finite numeric max_identity"):
        validate_split_payload(payload)
    payload = _valid_payload()
    payload["leak_check"]["strict_probe"]["test_vs_train"]["max_identity"] = 0.5
    with pytest.raises(SplitError, match="0.0 when there are no hits"):
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
    assert payload["mmseqs"]["executable"]
    assert set(payload["platform"]) == {"system", "machine", "python"}
    assert all(payload["platform"].values())
    assert payload["mmseqs"]["parameters"]["cluster_mode"] == 0
    assert payload["mmseqs"]["parameters"]["cluster_reassign"] is True
    assert payload["leak_check"]["hits"] == 0
    assert payload["schema_version"] == 3
    assert payload["leak_check"]["holdout_cross_check"]["validation_vs_test"] == {
        "queries_with_hits": 0, "max_identity": 0.0,
    }
    assert payload["leak_check"]["strict_probe"]["sensitivity"] == 7.5
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
        executable: str,
        query: Path,
        target: Path,
        work_dir: Path,
        label: str,
        *,
        sensitivity: float = 5.0,
    ) -> list[tuple[str, float]]:
        nonlocal calls
        calls += 1
        return [("P2", 0.9)] if calls == 1 else []

    monkeypatch.setattr("pfpred.splits._search_hits", fake_search)
    result = _repair_and_check_leakage(
        "mmseqs", clusters, assignments, sequences, tmp_path,
    )

    assert assignments == {"P1": "train", "P2": "train", "P2c": "validation", "P3": "test"}
    assert result["hits"] == 0
    assert result["moved_clusters_to_train"] == 1
    assert result["moved_proteins_to_train"] == 2


def test_leakage_repair_records_cross_holdout_and_strict_probe(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    clusters = {"P1": ("P1",), "P2": ("P2",), "P3": ("P3",)}
    assignments = {"P1": "train", "P2": "validation", "P3": "test"}
    sequences = {name: "ACDE" for name in clusters}
    seen: list[tuple[str, float]] = []

    def fake_search(
        executable: str,
        query: Path,
        target: Path,
        work_dir: Path,
        label: str,
        *,
        sensitivity: float = 5.0,
    ) -> list[tuple[str, float]]:
        seen.append((label, sensitivity))
        if label.startswith("cross_validation"):
            return [("P2", 0.91), ("P2", 0.55)]
        if label.startswith("strict_validation"):
            return [("P2", 0.99)]
        if label.startswith("strict_test"):
            return [("P3", 0.42), ("P3", 0.31)]
        return []

    monkeypatch.setattr("pfpred.splits._search_hits", fake_search)
    result = _repair_and_check_leakage(
        "mmseqs", clusters, assignments, sequences, tmp_path,
    )

    cross = result["holdout_cross_check"]
    assert cross["validation_vs_test"] == {"queries_with_hits": 1, "max_identity": 0.91}
    assert cross["test_vs_validation"] == {"queries_with_hits": 0, "max_identity": 0.0}
    strict = result["strict_probe"]
    assert strict["sensitivity"] == 7.5
    assert strict["validation_vs_train"] == {"queries_with_hits": 1, "max_identity": 0.99}
    assert strict["test_vs_train"] == {"queries_with_hits": 1, "max_identity": 0.42}
    assert {label for label, sens in seen if sens == 7.5} == {
        "strict_validation_1", "strict_test_1",
    }
    assert all(sens == 5.0 for label, sens in seen if label.startswith("cross_"))


def _fasta_headers(text: str) -> set[str]:
    return {line[1:] for line in text.splitlines() if line.startswith(">")}


def test_probes_read_final_post_repair_fastas(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    clusters = {"P1": ("P1",), "P2": ("P2", "P2b"), "P3": ("P3",), "P4": ("P4",)}
    assignments = {"P1": "train", "P2": "validation", "P3": "test", "P4": "validation"}
    sequences = {name: "ACDE" for members in clusters.values() for name in members}
    captured: dict[str, tuple[str, str]] = {}

    def fake_search(
        executable: str,
        query: Path,
        target: Path,
        work_dir: Path,
        label: str,
        *,
        sensitivity: float = 5.0,
    ) -> list[tuple[str, float]]:
        captured[label] = (
            query.read_text(encoding="utf-8"),
            target.read_text(encoding="utf-8"),
        )
        return [("P2", 0.9)] if label == "validation_1" else []

    monkeypatch.setattr("pfpred.splits._search_hits", fake_search)
    result = _repair_and_check_leakage(
        "mmseqs", clusters, assignments, sequences, tmp_path,
    )

    assert result["moved_clusters_to_train"] == 1
    assert result["final_protein_counts"] == {"train": 3, "validation": 1, "test": 1}
    strict_query, strict_target = captured["strict_validation_2"]
    assert _fasta_headers(strict_query) == {"P4"}
    assert {"P1", "P2", "P2b"} <= _fasta_headers(strict_target)
    assert "P3" not in _fasta_headers(strict_target)


def test_probes_reject_unknown_query_ids(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    clusters = {"P1": ("P1",), "P2": ("P2",), "P3": ("P3",)}
    assignments = {"P1": "train", "P2": "validation", "P3": "test"}
    sequences = {name: "ACDE" for name in clusters}

    def fake_search(
        executable: str,
        query: Path,
        target: Path,
        work_dir: Path,
        label: str,
        *,
        sensitivity: float = 5.0,
    ) -> list[tuple[str, float]]:
        return [("P9", 0.9)] if label.startswith("strict_validation") else []

    monkeypatch.setattr("pfpred.splits._search_hits", fake_search)
    with pytest.raises(SplitError, match="unknown query"):
        _repair_and_check_leakage(
            "mmseqs", clusters, assignments, sequences, tmp_path,
        )


def test_validate_split_payload_rejects_empty_partition() -> None:
    payload = _valid_payload()
    payload["assignments"].pop("P4")
    payload["split"]["protein_counts"]["test"] = 0
    payload["split"]["cluster_counts"]["test"] = 0
    payload["split"]["total_proteins"] = 3
    payload["split"]["total_clusters"] = 2
    with pytest.raises(SplitError, match="empty split"):
        validate_split_payload(payload)
