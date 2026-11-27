"""MMseqs2 clustering and leakage-safe protein-level dataset splits."""

from __future__ import annotations

import json
import random
import subprocess
import tempfile
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path

from pfpred.paths import sha256_file
from pfpred.uniprot import read_fasta

MIN_SEQ_ID = 0.30
MIN_COVERAGE = 0.80
COVERAGE_MODE = 0
SPLIT_RATIOS = {"train": 0.80, "validation": 0.10, "test": 0.10}
SPLIT_NAMES = tuple(SPLIT_RATIOS)


class SplitError(ValueError):
    """Raised when cluster output cannot produce a valid leakage-safe split."""


def read_cluster_tsv(
    path: Path, expected_accessions: set[str]
) -> dict[str, tuple[str, ...]]:
    """Read MMseqs2 representative/member TSV and validate exact input coverage.

    The expected set comes from the source FASTA. Every input accession must
    occur exactly once as a member, and every representative must be a member
    of its own cluster (including singleton clusters).
    """
    clusters: dict[str, list[str]] = {}
    seen_members: set[str] = set()

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 2 or not all(fields):
                raise SplitError(
                    f"{path}:{line_number}: expected representative<TAB>member, got {line!r}"
                )
            representative, member = fields
            if representative not in expected_accessions:
                raise SplitError(
                    f"{path}:{line_number}: unknown representative {representative!r}"
                )
            if member not in expected_accessions:
                raise SplitError(f"{path}:{line_number}: unknown member {member!r}")
            if member in seen_members:
                raise SplitError(f"{path}:{line_number}: duplicate cluster member {member!r}")
            seen_members.add(member)
            clusters.setdefault(representative, []).append(member)

    if not clusters:
        raise SplitError(f"{path}: MMseqs2 produced an empty cluster file")

    for representative, members in clusters.items():
        if representative not in members:
            raise SplitError(
                f"{path}: representative {representative!r} is absent from its own cluster"
            )

    missing = expected_accessions - seen_members
    if missing:
        examples = ", ".join(sorted(missing)[:5])
        raise SplitError(
            f"{path}: {len(missing)} FASTA proteins are missing from cluster output "
            f"(for example: {examples})"
        )

    unexpected = seen_members - expected_accessions
    if unexpected:  # Defensive: members are checked while parsing.
        raise SplitError(f"{path}: unexpected cluster members: {sorted(unexpected)[:5]}")

    return {
        representative: tuple(sorted(members))
        for representative, members in sorted(clusters.items())
    }


def assign_cluster_splits(
    clusters: Mapping[str, Sequence[str]], *, seed: int = 42
) -> dict[str, str]:
    """Assign whole clusters to splits, balancing protein counts deterministically.

    Clusters are placed largest-first. A seeded shuffle determines the order of
    equal-sized clusters; each cluster goes to the split furthest below its
    target fraction. Cluster boundaries are never broken.
    """
    if not clusters:
        raise SplitError("cannot split an empty cluster collection")

    seen: set[str] = set()
    cluster_items: list[tuple[str, tuple[str, ...]]] = []
    for representative, raw_members in clusters.items():
        members = tuple(raw_members)
        if not members:
            raise SplitError(f"cluster {representative!r} is empty")
        if representative not in members:
            raise SplitError(
                f"representative {representative!r} is absent from its own cluster"
            )
        duplicates = seen.intersection(members)
        if duplicates:
            raise SplitError(f"proteins occur in more than one cluster: {sorted(duplicates)[:5]}")
        seen.update(members)
        cluster_items.append((representative, members))

    total = len(seen)
    targets = {name: total * ratio for name, ratio in SPLIT_RATIOS.items()}
    counts = {name: 0 for name in SPLIT_NAMES}

    rng = random.Random(seed)
    rng.shuffle(cluster_items)
    cluster_items.sort(key=lambda item: len(item[1]), reverse=True)

    cluster_to_split: dict[str, str] = {}
    for representative, members in cluster_items:
        destination = min(
            SPLIT_NAMES,
            key=lambda name: counts[name] / targets[name],
        )
        cluster_to_split[representative] = destination
        counts[destination] += len(members)

    return cluster_to_split


def _write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        temporary_path.replace(path)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def run_clustering(
    fasta_path: Path,
    output_path: Path,
    *,
    mmseqs_executable: str = "mmseqs",
    seed: int = 42,
) -> dict[str, object]:
    """Run MMseqs2 easy-cluster and write validated split assignments as JSON."""
    accessions: list[str] = []
    seen_accessions: set[str] = set()
    for protein in read_fasta(fasta_path):
        if not protein.sequence:
            raise SplitError(f"FASTA contains an empty sequence for {protein.accession!r}")
        if protein.accession in seen_accessions:
            raise SplitError(f"FASTA contains duplicate accession {protein.accession!r}")
        accessions.append(protein.accession)
        seen_accessions.add(protein.accession)

    if not accessions:
        raise SplitError(f"FASTA contains no proteins: {fasta_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    version_result = subprocess.run(
        [mmseqs_executable, "version"],
        check=True,
        capture_output=True,
        text=True,
    )
    mmseqs_version = (version_result.stdout or version_result.stderr).strip()
    if not mmseqs_version:
        raise SplitError("MMseqs2 version command returned no version string")

    with tempfile.TemporaryDirectory(
        prefix="pfpred-mmseqs-", dir=output_path.parent
    ) as temporary_directory:
        work_dir = Path(temporary_directory)
        prefix = work_dir / "human"
        tmp_dir = work_dir / "tmp"
        tmp_dir.mkdir()
        command = [
            mmseqs_executable,
            "easy-cluster",
            str(fasta_path),
            str(prefix),
            str(tmp_dir),
            "--min-seq-id",
            f"{MIN_SEQ_ID:.2f}",
            "-c",
            f"{MIN_COVERAGE:.2f}",
            "--cov-mode",
            str(COVERAGE_MODE),
        ]
        subprocess.run(command, check=True)
        cluster_path = Path(f"{prefix}_cluster.tsv")
        clusters = read_cluster_tsv(cluster_path, seen_accessions)
        cluster_to_split = assign_cluster_splits(clusters, seed=seed)

    protein_assignments = {
        accession: {
            "cluster": representative,
            "split": cluster_to_split[representative],
        }
        for representative, members in clusters.items()
        for accession in members
    }
    protein_counts = Counter(
        assignment["split"] for assignment in protein_assignments.values()
    )
    cluster_counts = Counter(cluster_to_split.values())

    payload: dict[str, object] = {
        "schema_version": 1,
        "track": "A",
        "dataset": "Human Swiss-Prot 2026_04",
        "input_fasta": fasta_path.name,
        "input_sha256": sha256_file(fasta_path),
        "mmseqs": {
            "version": mmseqs_version,
            "workflow": "easy-cluster",
            "parameters": {
                "min_seq_id": MIN_SEQ_ID,
                "coverage": MIN_COVERAGE,
                "coverage_mode": COVERAGE_MODE,
            },
        },
        "split": {
            "seed": seed,
            "target_ratios": SPLIT_RATIOS,
            "protein_counts": {
                name: protein_counts.get(name, 0) for name in SPLIT_NAMES
            },
            "cluster_counts": {
                name: cluster_counts.get(name, 0) for name in SPLIT_NAMES
            },
            "total_proteins": len(protein_assignments),
            "total_clusters": len(clusters),
        },
        "assignments": dict(sorted(protein_assignments.items())),
    }
    _write_json(output_path, payload)
    return payload
