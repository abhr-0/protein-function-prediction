"""MMseqs2 clustering and leakage-checked protein dataset splits."""

from __future__ import annotations

import hashlib
import json
import os
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
CLUSTER_MODE = 0
CLUSTER_REASSIGN = True
SEARCH_SENSITIVITY = 5.0
SEARCH_MAX_SEQS = 1000
SPLIT_RATIOS = {"train": 0.80, "validation": 0.10, "test": 0.10}
SPLIT_NAMES = tuple(SPLIT_RATIOS)


class SplitError(ValueError):
    """Raised when clustering or split validation fails."""


def _is_windows_batch(executable: str) -> bool:
    return os.name == "nt" and Path(executable).suffix.lower() in {".bat", ".cmd"}


def _run_mmseqs(
    executable: str, arguments: Sequence[str], **kwargs: object
) -> subprocess.CompletedProcess[str]:
    """Run MMseqs2, using cmd.exe for Windows batch-file launchers."""
    command = (
        ["cmd.exe", "/d", "/c", "call", executable, *arguments]
        if _is_windows_batch(executable)
        else [executable, *arguments]
    )
    return subprocess.run(command, **kwargs)


def read_cluster_tsv(
    path: Path, expected_accessions: set[str]
) -> dict[str, tuple[str, ...]]:
    """Read representative/member TSV and require exact one-to-one FASTA coverage."""
    clusters: dict[str, list[str]] = {}
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
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
                raise SplitError(f"{path}:{line_number}: unknown representative {representative!r}")
            if member not in expected_accessions:
                raise SplitError(f"{path}:{line_number}: unknown member {member!r}")
            if member in seen:
                raise SplitError(f"{path}:{line_number}: duplicate cluster member {member!r}")
            seen.add(member)
            clusters.setdefault(representative, []).append(member)

    if not clusters:
        raise SplitError(f"{path}: MMseqs2 produced an empty cluster file")
    for representative, members in clusters.items():
        if representative not in members:
            raise SplitError(
                f"{path}: representative {representative!r} is absent from its own cluster"
            )
    missing = expected_accessions - seen
    if missing:
        raise SplitError(
            f"{path}: {len(missing)} FASTA proteins are missing from cluster output "
            f"(for example: {', '.join(sorted(missing)[:5])})"
        )
    return {rep: tuple(sorted(members)) for rep, members in sorted(clusters.items())}


def assign_cluster_splits(
    clusters: Mapping[str, Sequence[str]], *, seed: int = 42
) -> dict[str, str]:
    """Assign whole clusters deterministically, independent of mapping insertion order."""
    if not clusters:
        raise SplitError("cannot split an empty cluster collection")
    seen: set[str] = set()
    items: list[tuple[str, tuple[str, ...]]] = []
    for representative, raw_members in clusters.items():
        members = tuple(sorted(raw_members))
        if not members:
            raise SplitError(f"cluster {representative!r} is empty")
        if representative not in members:
            raise SplitError(f"representative {representative!r} is absent from its own cluster")
        overlap = seen.intersection(members)
        if overlap:
            raise SplitError(f"proteins occur in more than one cluster: {sorted(overlap)[:5]}")
        seen.update(members)
        items.append((representative, members))

    def tie_key(item: tuple[str, tuple[str, ...]]) -> str:
        return hashlib.sha256(f"{seed}:{item[0]}".encode()).hexdigest()

    items.sort(key=lambda item: (-len(item[1]), tie_key(item)))
    total = len(seen)
    targets = {name: total * ratio for name, ratio in SPLIT_RATIOS.items()}
    counts = dict.fromkeys(SPLIT_NAMES, 0)
    result: dict[str, str] = {}
    for representative, members in items:
        destination = min(SPLIT_NAMES, key=lambda name: (counts[name] / targets[name], name))
        result[representative] = destination
        counts[destination] += len(members)
    if any(counts[name] == 0 for name in SPLIT_NAMES):
        raise SplitError("cluster collection is too small to populate all three splits")
    return result


def _write_fasta(path: Path, accessions: Sequence[str], sequences: Mapping[str, str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for accession in sorted(accessions):
            handle.write(f">{accession}\n{sequences[accession]}\n")


def _search_hits(
    executable: str,
    query_path: Path,
    target_path: Path,
    work_dir: Path,
    label: str,
) -> set[str]:
    result_path = work_dir / f"{label}.tsv"
    tmp_path = work_dir / f"{label}_tmp"
    tmp_path.mkdir()
    _run_mmseqs(
        executable,
        [
            "easy-search", str(query_path), str(target_path), str(result_path), str(tmp_path),
            "--min-seq-id", f"{MIN_SEQ_ID:.2f}", "-c", f"{MIN_COVERAGE:.2f}",
            "--cov-mode", str(COVERAGE_MODE), "-s", str(SEARCH_SENSITIVITY),
            "--max-seqs", str(SEARCH_MAX_SEQS), "--format-output", "query,target", "-v", "1",
        ],
        check=True,
    )
    hits: set[str] = set()
    if result_path.exists():
        with result_path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 2 or not all(fields):
                    raise SplitError(f"{result_path}:{line_number}: malformed easy-search output")
                hits.add(fields[0])
    return hits


def _repair_and_check_leakage(
    executable: str,
    clusters: Mapping[str, Sequence[str]],
    cluster_to_split: dict[str, str],
    sequences: Mapping[str, str],
    work_dir: Path,
) -> dict[str, object]:
    """Move any query cluster with a train hit into train, repeating to a fixed point."""
    moved_clusters = 0
    moved_proteins = 0
    total_hits = 0
    rounds = 0
    while True:
        rounds += 1
        split_accessions = {
            name: sorted(
                member
                for representative, members in clusters.items()
                if cluster_to_split[representative] == name
                for member in members
            )
            for name in SPLIT_NAMES
        }
        if not split_accessions["train"]:
            raise SplitError("leakage repair left the training split empty")
        train_path = work_dir / "train.fasta"
        _write_fasta(train_path, split_accessions["train"], sequences)
        offending: set[str] = set()
        round_hits = 0
        for split_name in ("validation", "test"):
            query_accessions = split_accessions[split_name]
            if not query_accessions:
                raise SplitError(f"leakage repair emptied the {split_name} split")
            query_path = work_dir / f"{split_name}.fasta"
            _write_fasta(query_path, query_accessions, sequences)
            hits = _search_hits(
                executable, query_path, train_path, work_dir, f"{split_name}_{rounds}"
            )
            unknown = hits - set(query_accessions)
            if unknown:
                raise SplitError(
                    "MMseqs2 leakage search returned unknown query IDs: "
                    f"{sorted(unknown)[:5]}"
                )
            round_hits += len(hits)
            representative_by_member = {
                member: representative
                for representative, members in clusters.items()
                for member in members
            }
            offending.update(representative_by_member[accession] for accession in hits)
        total_hits += round_hits
        if not offending:
            if any(not split_accessions[name] for name in SPLIT_NAMES):
                raise SplitError("leakage check produced an empty split")
            return {
                "method": "MMseqs2 easy-search against current training split",
                "min_seq_id": MIN_SEQ_ID,
                "coverage": MIN_COVERAGE,
                "coverage_mode": COVERAGE_MODE,
                "sensitivity": SEARCH_SENSITIVITY,
                "max_seqs": SEARCH_MAX_SEQS,
                "hits": 0,
                "rounds": rounds,
                "hits_observed_before_repair": total_hits,
                "moved_clusters_to_train": moved_clusters,
                "moved_proteins_to_train": moved_proteins,
                "final_protein_counts": {name: len(split_accessions[name]) for name in SPLIT_NAMES},
            }
        before = len(offending)
        moved_clusters += before
        moved_proteins += sum(len(clusters[rep]) for rep in offending)
        for representative in offending:
            cluster_to_split[representative] = "train"
        if before == 0:
            raise SplitError("leakage repair made no progress")


def validate_split_payload(payload: Mapping[str, object]) -> None:
    """Validate split JSON structure, counts, labels, and cluster indivisibility."""
    assignments = payload.get("assignments")
    split = payload.get("split")
    if not isinstance(assignments, dict) or not assignments:
        raise SplitError("split payload assignments must be a non-empty object")
    if not isinstance(split, dict):
        raise SplitError("split payload is missing split metadata")
    protein_counts = Counter()
    cluster_labels: dict[str, str] = {}
    cluster_counts = Counter()
    for accession, assignment in assignments.items():
        if not isinstance(accession, str) or not isinstance(assignment, dict):
            raise SplitError("malformed protein assignment")
        label, cluster = assignment.get("split"), assignment.get("cluster")
        if label not in SPLIT_NAMES or not isinstance(cluster, str) or not cluster:
            raise SplitError(f"invalid assignment for {accession!r}")
        previous = cluster_labels.setdefault(cluster, label)
        if previous != label:
            raise SplitError(f"cluster {cluster!r} is divided across splits")
        protein_counts[label] += 1
    for _cluster, label in cluster_labels.items():
        cluster_counts[label] += 1
    if any(protein_counts[name] == 0 for name in SPLIT_NAMES):
        raise SplitError("split payload contains an empty split")
    if sum(protein_counts.values()) != split.get("total_proteins"):
        raise SplitError("total_proteins does not match assignments")
    if len(cluster_labels) != split.get("total_clusters"):
        raise SplitError("total_clusters does not match assignments")
    expected_proteins = split.get("protein_counts")
    expected_clusters = split.get("cluster_counts")
    if not isinstance(expected_proteins, dict) or not isinstance(expected_clusters, dict):
        raise SplitError("split count metadata is malformed")
    if dict(protein_counts) != {name: expected_proteins.get(name, 0) for name in SPLIT_NAMES}:
        raise SplitError("protein_counts do not match assignments")
    if dict(cluster_counts) != {name: expected_clusters.get(name, 0) for name in SPLIT_NAMES}:
        raise SplitError("cluster_counts do not match assignments")
    leak_check = payload.get("leak_check")
    if not isinstance(leak_check, dict) or leak_check.get("hits") != 0:
        raise SplitError("split payload does not contain a successful zero-hit leakage check")


def load_split(path: Path, *, fasta_path: Path | None = None) -> dict[str, object]:
    """Load and validate split JSON; optionally verify its source FASTA digest."""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SplitError(f"could not read split JSON {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise SplitError("split JSON root must be an object")
    validate_split_payload(payload)
    if fasta_path is not None and payload.get("input_sha256") != sha256_file(fasta_path):
        raise SplitError("split JSON input_sha256 does not match the supplied FASTA")
    return payload


def _write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        temporary.replace(path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def run_clustering(
    fasta_path: Path,
    output_path: Path,
    *,
    mmseqs_executable: str = "mmseqs",
    seed: int = 42,
) -> dict[str, object]:
    """Cluster proteins, repair cross-split hits, and save validated provenance."""
    sequences: dict[str, str] = {}
    for protein in read_fasta(fasta_path):
        if not protein.sequence:
            raise SplitError(f"FASTA contains an empty sequence for {protein.accession!r}")
        if protein.accession in sequences:
            raise SplitError(f"FASTA contains duplicate accession {protein.accession!r}")
        sequences[protein.accession] = protein.sequence
    if not sequences:
        raise SplitError(f"FASTA contains no proteins: {fasta_path}")

    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    version_result = _run_mmseqs(
        mmseqs_executable, ["version"], check=True, capture_output=True, text=True,
    )
    version = (version_result.stdout or version_result.stderr).strip()
    if not version or "GITDIR-NOTFOUND" in version:
        version = f"unreported (executable: {mmseqs_executable})"

    with tempfile.TemporaryDirectory(prefix="pfpred-mmseqs-", dir=output_path.parent) as tmp:
        work_dir = Path(tmp)
        normalized_fasta = work_dir / "input.fasta"
        _write_fasta(normalized_fasta, list(sequences), sequences)
        prefix = work_dir / "human"
        cluster_tmp = work_dir / "cluster_tmp"
        cluster_tmp.mkdir()
        _run_mmseqs(
            mmseqs_executable,
            [
                "easy-cluster", str(normalized_fasta), str(prefix), str(cluster_tmp),
                "--min-seq-id", f"{MIN_SEQ_ID:.2f}", "-c", f"{MIN_COVERAGE:.2f}",
                "--cov-mode", str(COVERAGE_MODE), "--cluster-mode", str(CLUSTER_MODE),
                "--cluster-reassign", "1" if CLUSTER_REASSIGN else "0", "-v", "1",
            ],
            check=True,
        )
        clusters = read_cluster_tsv(Path(f"{prefix}_cluster.tsv"), set(sequences))
        cluster_to_split = assign_cluster_splits(clusters, seed=seed)
        leak_check = _repair_and_check_leakage(
            mmseqs_executable, clusters, cluster_to_split, sequences, work_dir,
        )

    assignments = {
        member: {"cluster": representative, "split": cluster_to_split[representative]}
        for representative, members in clusters.items()
        for member in members
    }
    protein_counts = Counter(item["split"] for item in assignments.values())
    cluster_counts = Counter(cluster_to_split.values())
    payload: dict[str, object] = {
        "schema_version": 2,
        "track": "A",
        "dataset": "Human Swiss-Prot 2026_04",
        "input_fasta": fasta_path.name,
        "input_sha256": sha256_file(fasta_path),
        "mmseqs": {
            "version": version,
            "workflow": "easy-cluster",
            "parameters": {
                "min_seq_id": MIN_SEQ_ID,
                "coverage": MIN_COVERAGE,
                "coverage_mode": COVERAGE_MODE,
                "cluster_mode": CLUSTER_MODE,
                "cluster_reassign": CLUSTER_REASSIGN,
                "search_sensitivity": SEARCH_SENSITIVITY,
                "search_max_seqs": SEARCH_MAX_SEQS,
            },
        },
        "split": {
            "seed": seed,
            "target_ratios": SPLIT_RATIOS,
            "protein_counts": {name: protein_counts.get(name, 0) for name in SPLIT_NAMES},
            "cluster_counts": {name: cluster_counts.get(name, 0) for name in SPLIT_NAMES},
            "total_proteins": len(assignments),
            "total_clusters": len(clusters),
        },
        "leak_check": leak_check,
        "assignments": dict(sorted(assignments.items())),
    }
    validate_split_payload(payload)
    _write_json(output_path, payload)
    return payload
