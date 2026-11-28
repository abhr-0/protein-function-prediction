"""Train-only Track A GO vocabularies and label coverage reports."""

from __future__ import annotations

import hashlib
import json
import tempfile
from collections import Counter
from collections.abc import Iterable, Mapping
from pathlib import Path

from pfpred.labels import CURATED, filter_annotations
from pfpred.paths import sha256_file
from pfpred.records import Annotation, VALID_ASPECTS
from pfpred.splits import SPLIT_NAMES, load_split
from pfpred.uniprot import read_annotations

DEFAULT_MIN_COUNT = 5
VOCABULARY_SCHEMA_VERSION = 1
ASPECT_NAMES = {"F": "MF", "P": "BP", "C": "CC"}


class VocabularyError(ValueError):
    """Raised when vocabulary inputs or split assignments are inconsistent."""


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


def build_vocabulary(
    annotations: Iterable[Annotation],
    assignments: Mapping[str, str],
    *,
    min_count: int = DEFAULT_MIN_COUNT,
) -> tuple[dict[str, dict[str, int]], dict[str, dict[str, object]]]:
    """Fit per-aspect vocabularies on train proteins and report split coverage.

    Counts are distinct protein-term pairs. Validation and test labels are
    inspected only after the vocabulary is fixed and cannot affect selection.
    """
    if isinstance(min_count, bool) or min_count < 1:
        raise VocabularyError("min_count must be an integer >= 1")
    invalid = set(assignments.values()) - set(SPLIT_NAMES)
    if invalid:
        raise VocabularyError(f"unknown split names in assignments: {sorted(invalid)}")

    grouped: dict[str, dict[str, set[str]]] = {aspect: {} for aspect in VALID_ASPECTS}
    unassigned_proteins: set[str] = set()
    for row in filter_annotations(annotations, CURATED):
        if row.entry not in assignments:
            unassigned_proteins.add(row.entry)
            continue
        grouped[row.aspect].setdefault(row.entry, set()).add(row.term)
    vocabularies: dict[str, dict[str, int]] = {}
    coverage: dict[str, dict[str, object]] = {}
    for aspect in sorted(VALID_ASPECTS):
        protein_terms = grouped[aspect]
        split_proteins = {
            name: {
                protein: terms for protein, terms in protein_terms.items()
                if assignments.get(protein) == name
            }
            for name in SPLIT_NAMES
        }
        train_counts: Counter[str] = Counter(
            term for terms in split_proteins["train"].values() for term in terms
        )
        selected = {
            term: train_counts[term]
            for term in sorted(train_counts)
            if train_counts[term] >= min_count
        }
        vocabularies[aspect] = selected

        split_report: dict[str, object] = {}
        for split_name in SPLIT_NAMES:
            rows = split_proteins[split_name]
            total_pairs = sum(len(terms) for terms in rows.values())
            covered_pairs = sum(
                len(terms.intersection(selected)) for terms in rows.values()
            )
            proteins_with_labels = len(rows)
            proteins_with_positives = sum(
                bool(terms.intersection(selected)) for terms in rows.values()
            )
            split_report[split_name] = {
                "proteins_with_curated_annotations": proteins_with_labels,
                "proteins_with_positives": proteins_with_positives,
                "proteins_dropped_zero_positives": proteins_with_labels - proteins_with_positives,
                "annotation_pairs": total_pairs,
                "annotation_pairs_covered": covered_pairs,
                "annotation_pair_coverage": (
                    covered_pairs / total_pairs if total_pairs else 0.0
                ),
            }
        coverage[aspect] = {
            "aspect": ASPECT_NAMES[aspect],
            "min_count": min_count,
            "terms_kept": len(selected),
            "train_positive_proteins": sum(
                bool(terms.intersection(selected))
                for terms in split_proteins["train"].values()
            ),
            "splits": split_report,
        }
    coverage["_unassigned_curated_annotation_proteins"] = len(unassigned_proteins)
    return vocabularies, coverage


def build_vocabulary_files(
    annotation_path: Path,
    split_path: Path,
    vocabulary_path: Path,
    coverage_path: Path,
    *,
    min_count: int = DEFAULT_MIN_COUNT,
) -> tuple[dict[str, object], dict[str, object]]:
    """Read Track A inputs, fit vocabulary, and write provenance + coverage JSON."""
    split_payload = load_split(split_path)
    if split_payload.get("track") != "A":
        raise VocabularyError("Task B requires a Track A split artifact")
    assignments_raw = split_payload.get("assignments")
    if not isinstance(assignments_raw, dict):
        raise VocabularyError("split artifact assignments are malformed")
    assignments: dict[str, str] = {}
    for accession, item in assignments_raw.items():
        if not isinstance(item, dict) or not isinstance(item.get("split"), str):
            raise VocabularyError(f"malformed split assignment for {accession!r}")
        assignments[accession] = item["split"]

    vocabularies, coverage = build_vocabulary(
        read_annotations(annotation_path), assignments, min_count=min_count
    )
    split_digest = sha256_file(split_path)
    annotation_digest = sha256_file(annotation_path)
    vocab_payload: dict[str, object] = {
        "schema_version": VOCABULARY_SCHEMA_VERSION,
        "track": "A",
        "dataset": split_payload.get("dataset"),
        "evidence_policy": {
            "name": CURATED.name,
            "excluded": sorted(CURATED.excluded),
        },
        "selection": {
            "min_count": min_count,
            "count_unit": "distinct training proteins per GO term",
            "ancestor_propagation": False,
            "fit_split": "train",
        },
        "provenance": {
            "split_file": split_path.name,
            "split_sha256": split_digest,
            "annotations_file": annotation_path.name,
            "annotations_sha256": annotation_digest,
        },
        "aspects": {
            aspect: {
                "name": ASPECT_NAMES[aspect],
                "term_count": len(terms),
                "terms": terms,
            }
            for aspect, terms in vocabularies.items()
        },
    }
    coverage_payload: dict[str, object] = {
        "schema_version": 1,
        "track": "A",
        "min_count": min_count,
        "split_sha256": split_digest,
        "annotations_sha256": annotation_digest,
        "unassigned_curated_annotation_proteins": coverage.pop(
            "_unassigned_curated_annotation_proteins"
        ),
        "aspects": coverage,
    }
    _write_json(vocabulary_path, vocab_payload)
    _write_json(coverage_path, coverage_payload)
    return vocab_payload, coverage_payload
