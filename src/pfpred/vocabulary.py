"""Train-only Track A GO vocabularies and label coverage reports."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping
from numbers import Integral
from pathlib import Path

from pfpred.io import write_json_atomic
from pfpred.labels import CURATED, filter_annotations, group_by_aspect
from pfpred.paths import sha256_file
from pfpred.records import VALID_ASPECTS, Annotation, DataFormatError
from pfpred.splits import SPLIT_NAMES, load_split
from pfpred.uniprot import read_annotations

DEFAULT_MIN_COUNT = 5
VOCABULARY_SCHEMA_VERSION = 1
COVERAGE_SCHEMA_VERSION = 1
ASPECT_NAMES = {"F": "MF", "P": "BP", "C": "CC"}


class VocabularyError(ValueError):
    """Raised when vocabulary inputs or split assignments are inconsistent."""


def build_vocabulary(
    annotations: Iterable[Annotation],
    assignments: Mapping[str, str],
    *,
    min_count: int = DEFAULT_MIN_COUNT,
) -> tuple[dict[str, dict[str, int]], dict[str, dict[str, object]], int]:
    """Fit per-aspect vocabularies on train proteins and report split coverage.

    Counts are distinct protein-term pairs. Validation and test labels are
    inspected only after the vocabulary is fixed and cannot affect selection.
    """
    if not isinstance(min_count, Integral) or isinstance(min_count, bool) or min_count < 1:
        raise VocabularyError("min_count must be an integer >= 1")
    min_count = int(min_count)
    invalid = set(assignments.values()) - set(SPLIT_NAMES)
    if invalid:
        raise VocabularyError(f"unknown split names in assignments: {sorted(invalid)}")

    curated = list(filter_annotations(annotations, CURATED))
    unassigned = {row.entry for row in curated if row.entry not in assignments}
    if unassigned:
        examples = ", ".join(sorted(unassigned)[:5])
        raise VocabularyError(
            f"{len(unassigned)} curated annotation proteins are absent from the split "
            f"(for example: {examples})"
        )

    try:
        grouped = group_by_aspect(curated)
    except DataFormatError as exc:
        raise VocabularyError(str(exc)) from exc

    vocabularies: dict[str, dict[str, int]] = {}
    coverage: dict[str, dict[str, object]] = {}
    for aspect in sorted(VALID_ASPECTS):
        protein_terms = grouped.get(aspect, {})
        split_proteins = {
            name: {
                protein: set(terms) for protein, terms in protein_terms.items()
                if assignments[protein] == name
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
            proteins_in_split = sum(
                label == split_name for label in assignments.values()
            )
            proteins_with_labels = len(rows)
            proteins_with_positives = sum(
                bool(terms.intersection(selected)) for terms in rows.values()
            )
            split_report[split_name] = {
                "proteins_in_split": proteins_in_split,
                "proteins_with_curated_annotations": proteins_with_labels,
                "proteins_without_aspect_label": proteins_in_split - proteins_with_labels,
                "proteins_dropped_by_min_count": proteins_with_labels - proteins_with_positives,
                "proteins_with_positives": proteins_with_positives,
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
    return vocabularies, coverage, len(unassigned)


def _validate_distinct_paths(
    annotation_path: Path,
    split_path: Path,
    vocabulary_path: Path,
    coverage_path: Path,
    fasta_path: Path,
) -> None:
    paths = {
        "annotations": annotation_path.expanduser().resolve(),
        "split": split_path.expanduser().resolve(),
        "fasta": fasta_path.expanduser().resolve(),
        "vocabulary": vocabulary_path.expanduser().resolve(),
        "coverage": coverage_path.expanduser().resolve(),
    }
    reverse: dict[Path, str] = {}
    for label, path in paths.items():
        if path in reverse:
            raise VocabularyError(
                f"path collision: {label} and {reverse[path]} both resolve to {path}"
            )
        reverse[path] = label


def build_vocabulary_files(
    annotation_path: Path,
    split_path: Path,
    vocabulary_path: Path,
    coverage_path: Path,
    *,
    fasta_path: Path,
    min_count: int = DEFAULT_MIN_COUNT,
) -> tuple[dict[str, object], dict[str, object]]:
    """Read Track A inputs, fit vocabulary, and write provenance + coverage JSON."""
    _validate_distinct_paths(
        annotation_path, split_path, vocabulary_path, coverage_path, fasta_path
    )
    split_payload = load_split(split_path, fasta_path=fasta_path)
    if split_payload.get("track") != "A":
        raise VocabularyError(
            f"vocabulary requires a Track A split artifact, got {split_payload.get('track')!r}"
        )
    assignments_raw = split_payload["assignments"]
    assignments = {
        accession: item["split"]
        for accession, item in assignments_raw.items()
    }

    vocabularies, coverage, unassigned_count = build_vocabulary(
        read_annotations(annotation_path), assignments, min_count=min_count
    )
    empty_aspects = [
        ASPECT_NAMES[aspect] for aspect, terms in vocabularies.items() if not terms
    ]
    if empty_aspects:
        raise VocabularyError(
            f"vocabulary is empty for aspect(s): {', '.join(empty_aspects)}"
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
            "min_count": int(min_count),
            "count_unit": "distinct training proteins per GO term",
            "ancestor_propagation": False,
            "fit_split": "train",
        },
        "provenance": {
            "split_file": split_path.name,
            "split_sha256": split_digest,
            "input_fasta": split_payload.get("input_fasta"),
            "input_fasta_sha256": split_payload.get("input_sha256"),
            "platform": split_payload.get("platform"),
            "mmseqs": split_payload.get("mmseqs"),
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
        "schema_version": COVERAGE_SCHEMA_VERSION,
        "track": "A",
        "min_count": int(min_count),
        "split_sha256": split_digest,
        "annotations_sha256": annotation_digest,
        "unassigned_curated_annotation_proteins": unassigned_count,
        "aspects": coverage,
    }
    write_json_atomic(vocabulary_path, vocab_payload)
    write_json_atomic(coverage_path, coverage_payload)
    return vocab_payload, coverage_payload
