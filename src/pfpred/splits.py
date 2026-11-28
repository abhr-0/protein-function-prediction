"""MMseqs2 clustering and leakage-checked protein dataset splits."""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from collections.abc import Collection, Mapping, Sequence
from pathlib import Path

from pfpred.io import write_json_atomic
from pfpred.paths import sha256_file
from pfpred.uniprot import read_fasta

MIN_SEQ_ID = 0.30
MIN_COVERAGE = 0.80
COVERAGE_MODE = 0
CLUSTER_MODE = 0
CLUSTER_REASSIGN = True
SEARCH_SENSITIVITY = 5.0
STRICT_SEARCH_SENSITIVITY = 7.5
SEARCH_MAX_SEQS = 1000
SPLIT_SCHEMA_VERSION = 3
SPLIT_RATIOS = {"train": 0.80, "validation": 0.10, "test": 0.10}
SPLIT_NAMES = tuple(SPLIT_RATIOS)
CROSS_CHECK_DIRECTIONS = ("validation_vs_test", "test_vs_validation")
STRICT_PROBE_DIRECTIONS = ("validation_vs_train", "test_vs_train")
DIRECTION_QUERY_SPLITS = {
    "validation_vs_test": "validation",
    "test_vs_validation": "test",
    "validation_vs_train": "validation",
    "test_vs_train": "test",
}
HEX_DIGITS = frozenset("0123456789abcdef")
SHA256_HEX_LENGTH = 64


class SplitError(ValueError):
    """Raised when clustering or split validation fails."""


def _is_windows_batch(executable: str) -> bool:
    return os.name == "nt" and Path(executable).suffix.lower() in {".bat", ".cmd"}


def _canonical_executable(executable: str) -> str:
    """Record the resolved path for an executable found on PATH."""
    found = shutil.which(executable)
    return str(Path(found).resolve()) if found else executable


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
        # SPLIT_NAMES order is the tie-break priority: train, then validation, then test.
        destination = min(SPLIT_NAMES, key=lambda name: counts[name] / targets[name])
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
    *,
    sensitivity: float = SEARCH_SENSITIVITY,
) -> list[tuple[str, float]]:
    """Run easy-search and return (query accession, fractional identity) rows."""
    result_path = work_dir / f"{label}.tsv"
    tmp_path = work_dir / f"{label}_tmp"
    try:
        tmp_path.mkdir()
    except FileExistsError as exc:
        raise SplitError(f"duplicate MMseqs2 search label {label!r}") from exc
    _run_mmseqs(
        executable,
        [
            "easy-search", str(query_path), str(target_path), str(result_path), str(tmp_path),
            "--min-seq-id", f"{MIN_SEQ_ID:.2f}", "-c", f"{MIN_COVERAGE:.2f}",
            "--cov-mode", str(COVERAGE_MODE), "-s", str(sensitivity),
            "--max-seqs", str(SEARCH_MAX_SEQS),
            "--format-output", "query,target,fident", "-v", "1",
        ],
        check=True,
    )
    hits: list[tuple[str, float]] = []
    if result_path.exists():
        with result_path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 3 or not all(fields):
                    raise SplitError(f"{result_path}:{line_number}: malformed easy-search output")
                try:
                    identity = float(fields[2])
                except ValueError as exc:
                    raise SplitError(
                        f"{result_path}:{line_number}: non-numeric identity {fields[2]!r}"
                    ) from exc
                if not math.isfinite(identity):
                    raise SplitError(
                        f"{result_path}:{line_number}: non-finite identity {fields[2]!r}"
                    )
                hits.append((fields[0], identity))
    return hits


def _summarize_hits(
    hits: Sequence[tuple[str, float]], expected_queries: Collection[str]
) -> dict[str, object]:
    """Summarize hits as distinct query proteins with a hit plus the best identity."""
    unknown = {query for query, _identity in hits} - set(expected_queries)
    if unknown:
        raise SplitError(f"MMseqs2 search returned unknown query IDs: {sorted(unknown)[:5]}")
    return {
        "queries_with_hits": len({query for query, _identity in hits}),
        "max_identity": max((identity for _query, identity in hits), default=0.0),
    }


def _probe_pair(
    executable: str,
    query_path: Path,
    target_path: Path,
    work_dir: Path,
    label: str,
    *,
    sensitivity: float,
    expected_queries: Collection[str],
) -> dict[str, object]:
    hits = _search_hits(
        executable, query_path, target_path, work_dir, label, sensitivity=sensitivity
    )
    return _summarize_hits(hits, expected_queries)


def _holdout_cross_check(
    executable: str, work_dir: Path, split_accessions: Mapping[str, Sequence[str]]
) -> dict[str, object]:
    """Search the two holdouts against each other; recorded evidence, never repaired."""
    return {
        "method": "MMseqs2 easy-search validation against test and test against validation",
        "min_seq_id": MIN_SEQ_ID,
        "coverage": MIN_COVERAGE,
        "coverage_mode": COVERAGE_MODE,
        "sensitivity": SEARCH_SENSITIVITY,
        "max_seqs": SEARCH_MAX_SEQS,
        "repaired": False,
        "validation_vs_test": _probe_pair(
            executable, work_dir / "validation.fasta", work_dir / "test.fasta",
            work_dir, "cross_validation_vs_test", sensitivity=SEARCH_SENSITIVITY,
            expected_queries=split_accessions["validation"],
        ),
        "test_vs_validation": _probe_pair(
            executable, work_dir / "test.fasta", work_dir / "validation.fasta",
            work_dir, "cross_test_vs_validation", sensitivity=SEARCH_SENSITIVITY,
            expected_queries=split_accessions["test"],
        ),
    }


def _strict_probe(
    executable: str,
    work_dir: Path,
    round_label: str,
    split_accessions: Mapping[str, Sequence[str]],
) -> dict[str, object]:
    """Record residual holdout-to-train hits at a stricter sensitivity than the repair."""
    return {
        "method": "MMseqs2 easy-search holdouts against train at stricter sensitivity",
        "min_seq_id": MIN_SEQ_ID,
        "coverage": MIN_COVERAGE,
        "coverage_mode": COVERAGE_MODE,
        "sensitivity": STRICT_SEARCH_SENSITIVITY,
        "max_seqs": SEARCH_MAX_SEQS,
        "repaired": False,
        "validation_vs_train": _probe_pair(
            executable, work_dir / "validation.fasta", work_dir / "train.fasta",
            work_dir, f"strict_validation_{round_label}",
            sensitivity=STRICT_SEARCH_SENSITIVITY,
            expected_queries=split_accessions["validation"],
        ),
        "test_vs_train": _probe_pair(
            executable, work_dir / "test.fasta", work_dir / "train.fasta",
            work_dir, f"strict_test_{round_label}", sensitivity=STRICT_SEARCH_SENSITIVITY,
            expected_queries=split_accessions["test"],
        ),
    }


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
    representative_by_member = {
        member: representative
        for representative, members in clusters.items()
        for member in members
    }
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
            queries = {query for query, _identity in hits}
            unknown = queries - set(query_accessions)
            if unknown:
                raise SplitError(
                    "MMseqs2 leakage search returned unknown query IDs: "
                    f"{sorted(unknown)[:5]}"
                )
            round_hits += len(queries)
            offending.update(representative_by_member[accession] for accession in queries)
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
                "holdout_cross_check": _holdout_cross_check(
                    executable, work_dir, split_accessions
                ),
                "strict_probe": _strict_probe(
                    executable, work_dir, str(rounds), split_accessions
                ),
            }
        before = len(offending)
        moved_clusters += before
        moved_proteins += sum(len(clusters[rep]) for rep in offending)
        for representative in offending:
            cluster_to_split[representative] = "train"


def _is_number(value: object) -> bool:
    """Return True for real numbers, excluding booleans."""
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _require_number(section: Mapping[str, object], key: str, context: str) -> float:
    """Return a finite numeric field, rejecting booleans and strings."""
    value = section.get(key)
    if not _is_number(value) or not math.isfinite(value):
        raise SplitError(f"{context}.{key} must be a finite number")
    return float(value)


def _require_parameter(
    section: Mapping[str, object], key: str, expected: float, context: str
) -> None:
    """Require a recorded parameter to equal the constant used by the repair gate."""
    value = _require_number(section, key, context)
    if value != expected:
        raise SplitError(f"{context}.{key} must be {expected!r}, got {value!r}")


def _require_integer(
    section: Mapping[str, object], key: str, context: str, *, minimum: int
) -> int:
    """Return an integer field, rejecting booleans and values below the minimum."""
    value = section.get(key)
    if not isinstance(value, int) or isinstance(value, bool):
        raise SplitError(f"{context}.{key} must be an integer")
    if value < minimum:
        raise SplitError(f"{context}.{key} must be at least {minimum}")
    return value


def _require_count(section: Mapping[str, object], key: str, context: str) -> int:
    """Return a non-negative integer count."""
    return _require_integer(section, key, context, minimum=0)


def _require_bool(
    section: Mapping[str, object], key: str, expected: bool, context: str
) -> None:
    if section.get(key) is not expected:
        raise SplitError(f"{context}.{key} must be {str(expected).lower()}")


def _require_gate_parameters(section: Mapping[str, object], context: str) -> None:
    """Require the search thresholds and max-seqs to match the repair gate."""
    _require_parameter(section, "min_seq_id", MIN_SEQ_ID, context)
    _require_parameter(section, "coverage", MIN_COVERAGE, context)
    _require_parameter(section, "coverage_mode", COVERAGE_MODE, context)
    _require_parameter(section, "max_seqs", SEARCH_MAX_SEQS, context)


def _require_cluster_parameters(mmseqs: Mapping[str, object], context: str) -> None:
    """Require the recorded easy-cluster parameters to match the module constants."""
    parameters = mmseqs.get("parameters")
    if not isinstance(parameters, dict):
        raise SplitError(f"{context}.parameters must be an object")
    context = f"{context}.parameters"
    _require_parameter(parameters, "min_seq_id", MIN_SEQ_ID, context)
    _require_parameter(parameters, "coverage", MIN_COVERAGE, context)
    _require_parameter(parameters, "coverage_mode", COVERAGE_MODE, context)
    _require_parameter(parameters, "cluster_mode", CLUSTER_MODE, context)
    _require_parameter(parameters, "search_sensitivity", SEARCH_SENSITIVITY, context)
    _require_parameter(parameters, "search_max_seqs", SEARCH_MAX_SEQS, context)
    _require_bool(parameters, "cluster_reassign", CLUSTER_REASSIGN, context)


def _require_false(section: Mapping[str, object], key: str, context: str) -> None:
    if section.get(key) is not False:
        raise SplitError(f"{context}.{key} must be false")


def _validate_summary(
    section: Mapping[str, object],
    direction: str,
    name: str,
    protein_counts: Mapping[str, int],
) -> None:
    context = f"leak_check.{name}.{direction}"
    summary = section.get(direction)
    if not isinstance(summary, dict):
        raise SplitError(f"leak_check.{name} is missing the {direction!r} summary")
    queries_with_hits = summary.get("queries_with_hits")
    if not isinstance(queries_with_hits, int) or isinstance(queries_with_hits, bool):
        raise SplitError(f"{context} needs an integer queries_with_hits")
    if queries_with_hits < 0:
        raise SplitError(f"{context} needs a non-negative queries_with_hits")
    query_split = DIRECTION_QUERY_SPLITS[direction]
    if queries_with_hits > protein_counts[query_split]:
        raise SplitError(
            f"{context} queries_with_hits exceeds the {query_split} protein count"
        )
    max_identity = summary.get("max_identity")
    if not _is_number(max_identity) or not math.isfinite(max_identity):
        raise SplitError(f"{context} needs a finite numeric max_identity")
    if not 0.0 <= max_identity <= 1.0:
        raise SplitError(f"{context} needs max_identity within [0.0, 1.0]")
    if queries_with_hits == 0 and max_identity != 0.0:
        raise SplitError(f"{context} needs max_identity 0.0 when there are no hits")


def validate_split_payload(payload: Mapping[str, object]) -> None:
    """Validate split JSON structure, counts, labels, and cluster indivisibility."""
    if payload.get("schema_version") != SPLIT_SCHEMA_VERSION:
        raise SplitError(f"split payload schema_version must be {SPLIT_SCHEMA_VERSION}")
    input_sha256 = payload.get("input_sha256")
    if (
        not isinstance(input_sha256, str)
        or len(input_sha256) != SHA256_HEX_LENGTH
        or not set(input_sha256) <= HEX_DIGITS
    ):
        raise SplitError("split payload input_sha256 must be a 64-character lowercase hex digest")
    platform_info = payload.get("platform")
    if not isinstance(platform_info, dict) or not all(
        isinstance(platform_info.get(key), str) and platform_info.get(key)
        for key in ("system", "machine", "python")
    ):
        raise SplitError("split payload platform needs non-empty system, machine, and python")
    mmseqs = payload.get("mmseqs")
    if not isinstance(mmseqs, dict):
        raise SplitError("split payload is missing mmseqs provenance")
    executable = mmseqs.get("executable")
    if not isinstance(executable, str) or not executable:
        raise SplitError("split payload mmseqs.executable must be a non-empty string")
    version = mmseqs.get("version")
    if not isinstance(version, str) or not version.strip():
        raise SplitError("split payload mmseqs.version must be a non-empty string")
    if mmseqs.get("workflow") != "easy-cluster":
        raise SplitError("split payload mmseqs.workflow must be 'easy-cluster'")
    _require_cluster_parameters(mmseqs, "split payload mmseqs")
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
    if not isinstance(leak_check, dict):
        raise SplitError("split payload is missing the leak_check section")
    hits = leak_check.get("hits")
    if not isinstance(hits, int) or isinstance(hits, bool) or hits != 0:
        raise SplitError("split payload does not contain a successful zero-hit leakage check")
    _require_parameter(leak_check, "sensitivity", SEARCH_SENSITIVITY, "leak_check")
    _require_gate_parameters(leak_check, "leak_check")
    _require_integer(leak_check, "rounds", "leak_check", minimum=1)
    for key in (
        "hits_observed_before_repair",
        "moved_clusters_to_train",
        "moved_proteins_to_train",
    ):
        _require_count(leak_check, key, "leak_check")
    final_counts = leak_check.get("final_protein_counts")
    if not isinstance(final_counts, dict) or set(final_counts) != set(SPLIT_NAMES):
        raise SplitError(
            "leak_check.final_protein_counts must contain exactly "
            f"{', '.join(SPLIT_NAMES)}"
        )
    for name in SPLIT_NAMES:
        _require_count(final_counts, name, "leak_check.final_protein_counts")
    if final_counts != expected_proteins:
        raise SplitError(
            "leak_check.final_protein_counts do not match split.protein_counts"
        )
    validated_counts = {name: protein_counts[name] for name in SPLIT_NAMES}
    cross_check = leak_check.get("holdout_cross_check")
    if not isinstance(cross_check, dict):
        raise SplitError("leak_check is missing the holdout_cross_check section")
    _require_parameter(
        cross_check, "sensitivity", SEARCH_SENSITIVITY, "leak_check.holdout_cross_check"
    )
    _require_gate_parameters(cross_check, "leak_check.holdout_cross_check")
    _require_false(cross_check, "repaired", "leak_check.holdout_cross_check")
    for direction in CROSS_CHECK_DIRECTIONS:
        _validate_summary(cross_check, direction, "holdout_cross_check", validated_counts)
    strict_probe = leak_check.get("strict_probe")
    if not isinstance(strict_probe, dict):
        raise SplitError("leak_check is missing the strict_probe section")
    if _require_number(
        strict_probe, "sensitivity", "leak_check.strict_probe"
    ) <= SEARCH_SENSITIVITY:
        raise SplitError(
            "leak_check.strict_probe.sensitivity must exceed the repair sensitivity "
            f"{SEARCH_SENSITIVITY!r}"
        )
    _require_gate_parameters(strict_probe, "leak_check.strict_probe")
    _require_false(strict_probe, "repaired", "leak_check.strict_probe")
    for direction in STRICT_PROBE_DIRECTIONS:
        _validate_summary(strict_probe, direction, "strict_probe", validated_counts)


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


_write_json = write_json_atomic

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
        prefix = work_dir / "clusters"
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
        "schema_version": SPLIT_SCHEMA_VERSION,
        "track": "A",
        "dataset": "Human Swiss-Prot 2026_04",
        "input_fasta": fasta_path.name,
        "input_sha256": sha256_file(fasta_path),
        "platform": {
            "system": sys.platform,
            "machine": platform.machine(),
            "python": platform.python_version(),
        },
        "mmseqs": {
            "version": version,
            "executable": _canonical_executable(mmseqs_executable),
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
