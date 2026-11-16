"""Target GO vocabulary selection over per-protein term collections.

Selecting the target vocabulary *after* the train/test split is essential:
using term frequencies from the full dataset leaks test-set information into
the label space.
"""

from collections import Counter
from collections.abc import Collection, Iterable, Sequence


def count_protein_terms(protein_terms: Iterable[Collection[str]]) -> Counter[str]:
    """Count how many proteins are annotated with each GO term.

    Duplicate terms within one protein are counted once.
    """
    counts: Counter[str] = Counter()
    for terms in protein_terms:
        counts.update(set(terms))
    return counts


def select_target_terms(
    protein_terms: Iterable[Collection[str]],
    *,
    min_count: int = 1,
    top_k: int | None = None,
) -> list[str]:
    """Select target terms by protein frequency.

    Args:
        protein_terms: One collection of GO terms per protein.
        min_count: Minimum number of proteins a term must appear in.
        top_k: Optional cap on the vocabulary size.

    Returns:
        Term IDs ordered by descending frequency, then lexicographically for
        deterministic tie-breaking.
    """
    if min_count < 1:
        raise ValueError(f"min_count must be >= 1, got {min_count}")
    if top_k is not None and top_k < 1:
        raise ValueError(f"top_k must be >= 1 or None, got {top_k}")

    counts = count_protein_terms(protein_terms)
    ranked = sorted(
        ((term, count) for term, count in counts.items() if count >= min_count),
        key=lambda item: (-item[1], item[0]),
    )
    if top_k is not None:
        ranked = ranked[:top_k]
    return [term for term, _ in ranked]


def filter_terms(
    protein_terms: Sequence[Collection[str]],
    vocabulary: Collection[str],
) -> list[list[str]]:
    """Restrict each protein's terms to the target vocabulary.

    Returns sorted, de-duplicated term lists so downstream target matrices are
    deterministic.
    """
    vocabulary_set = set(vocabulary)
    return [
        [term for term in sorted(set(terms)) if term in vocabulary_set]
        for terms in protein_terms
    ]


def vocabulary_coverage(
    protein_terms: Sequence[Collection[str]],
    vocabulary: Collection[str],
) -> dict[str, float | int]:
    """Report how much of the observed label space the vocabulary covers.

    ``annotation_coverage`` is the fraction of all observed protein-term pairs
    represented by the vocabulary. It caps achievable recall when the
    vocabulary is truncated.
    """
    vocabulary_set = set(vocabulary)
    total_pairs = 0
    covered_pairs = 0
    covered_proteins = 0

    for terms in protein_terms:
        unique_terms = set(terms)
        total_pairs += len(unique_terms)
        matched = len(unique_terms & vocabulary_set)
        covered_pairs += matched
        if matched:
            covered_proteins += 1

    annotation_coverage = covered_pairs / total_pairs if total_pairs else 0.0
    protein_coverage = covered_proteins / len(protein_terms) if protein_terms else 0.0
    return {
        "proteins": len(protein_terms),
        "proteins_with_target": covered_proteins,
        "protein_coverage": round(protein_coverage, 6),
        "annotation_pairs": total_pairs,
        "covered_pairs": covered_pairs,
        "annotation_coverage": round(annotation_coverage, 6),
    }
