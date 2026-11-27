from __future__ import annotations

from pathlib import Path

import pytest

from pfpred.splits import (
    SplitError,
    assign_cluster_splits,
    read_cluster_tsv,
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
    assert all(split in {"train", "validation", "test"} for split in first.values())


def test_assign_cluster_splits_rejects_overlapping_clusters() -> None:
    clusters = {"P1": ("P1", "P2"), "P2": ("P2", "P3")}

    with pytest.raises(SplitError, match="more than one cluster"):
        assign_cluster_splits(clusters)


def test_assign_cluster_splits_rejects_empty_input() -> None:
    with pytest.raises(SplitError, match="empty cluster collection"):
        assign_cluster_splits({})
