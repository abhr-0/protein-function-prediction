"""Tests for train-only target vocabulary selection."""

import pytest

from src.data.targets import (
    count_protein_terms,
    filter_terms,
    select_target_terms,
    vocabulary_coverage,
)


def test_count_protein_terms_deduplicates_within_protein():
    """A term repeated in one protein should count as one protein."""
    counts = count_protein_terms(
        [
            ["GO:1", "GO:1", "GO:2"],
            ["GO:1"],
            [],
        ]
    )
    assert counts == {"GO:1": 2, "GO:2": 1}


def test_select_target_terms_orders_by_frequency_then_id():
    """Ties must be deterministic for reproducible experiments."""
    terms = [["GO:B", "GO:A"], ["GO:A"], ["GO:B"], ["GO:C"]]
    assert select_target_terms(terms) == ["GO:A", "GO:B", "GO:C"]


def test_select_target_terms_applies_min_count_and_top_k():
    """Frequency thresholds and the vocabulary cap should both apply."""
    terms = [["GO:A", "GO:B"], ["GO:A", "GO:C"], ["GO:A"], ["GO:D"]]
    selected = select_target_terms(terms, min_count=2, top_k=1)
    assert selected == ["GO:A"]


def test_select_target_terms_rejects_invalid_bounds():
    """Invalid thresholds should fail loudly rather than silently change data."""
    with pytest.raises(ValueError):
        select_target_terms([["GO:1"]], min_count=0)
    with pytest.raises(ValueError):
        select_target_terms([["GO:1"]], top_k=0)


def test_filter_terms_restricts_sorts_and_deduplicates():
    """Filtering should produce deterministic lists restricted to the vocabulary."""
    filtered = filter_terms(
        [["GO:2", "GO:1", "GO:2", "GO:9"], ["GO:9"]],
        vocabulary=["GO:1", "GO:2"],
    )
    assert filtered == [["GO:1", "GO:2"], []]


def test_vocabulary_coverage_reports_pair_and_protein_coverage():
    """Coverage should expose how much recall the vocabulary caps."""
    coverage = vocabulary_coverage(
        [["GO:1", "GO:2"], ["GO:2", "GO:3"], ["GO:9"]],
        vocabulary=["GO:1", "GO:2"],
    )
    assert coverage["proteins"] == 3
    assert coverage["proteins_with_target"] == 2
    assert coverage["protein_coverage"] == 0.666667
    assert coverage["annotation_pairs"] == 5
    assert coverage["covered_pairs"] == 3
    assert coverage["annotation_coverage"] == 0.6


def test_vocabulary_coverage_handles_empty_input():
    """Empty data should not divide by zero."""
    coverage = vocabulary_coverage([], vocabulary=["GO:1"])
    assert coverage["proteins"] == 0
    assert coverage["annotation_coverage"] == 0.0
