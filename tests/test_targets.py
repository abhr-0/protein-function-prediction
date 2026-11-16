"""Tests for train-only target vocabulary selection."""

import pandas as pd
import pytest

from src.data.targets import (
    apply_vocabulary,
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


def _frame(terms_per_protein: list[list[str]]) -> pd.DataFrame:
    """Build a minimal dataset frame shaped like build_annotated_dataset output."""
    return pd.DataFrame(
        {
            "EntryID": [f"P{index}" for index in range(len(terms_per_protein))],
            "terms": terms_per_protein,
        }
    )


def test_apply_vocabulary_drops_zero_positive_proteins():
    """Proteins with no in-vocabulary term must be removed, with stats."""
    frame = _frame([["GO:1", "GO:9"], ["GO:9"], ["GO:2"]])
    kept, stats = apply_vocabulary(frame, ["GO:1", "GO:2"])

    assert list(kept["EntryID"]) == ["P0", "P2"]
    assert kept["terms"].tolist() == [["GO:1"], ["GO:2"]]
    assert stats["proteins_in"] == 3
    assert stats["proteins_kept"] == 2
    assert stats["proteins_dropped"] == 1
    assert stats["annotation_pairs_in"] == 4
    assert stats["annotation_pairs_kept"] == 2
    assert stats["min_terms"] == 1


def test_apply_vocabulary_respects_min_terms():
    """min_terms > 1 should drop proteins with too few in-vocabulary terms."""
    frame = _frame([["GO:1", "GO:2"], ["GO:1"], ["GO:3"]])

    kept_one, _ = apply_vocabulary(frame, ["GO:1", "GO:2"], min_terms=1)
    kept_two, stats = apply_vocabulary(frame, ["GO:1", "GO:2"], min_terms=2)

    assert list(kept_one["EntryID"]) == ["P0", "P1"]
    assert list(kept_two["EntryID"]) == ["P0"]
    assert stats["proteins_dropped"] == 2


def test_apply_vocabulary_rejects_invalid_min_terms():
    """A zero or negative threshold would silently keep everything."""
    with pytest.raises(ValueError):
        apply_vocabulary(_frame([["GO:1"]]), ["GO:1"], min_terms=0)


def test_apply_vocabulary_reports_test_proteins_lost_to_vocabulary():
    """Test proteins outside the train vocabulary are dropped, and counted.

    This is the documented consequence of choosing the vocabulary on the train
    split only: those proteins fall outside the evaluated label space, so the
    drop count must stay visible rather than looking like full coverage.
    """
    train = _frame([["GO:1"], ["GO:1", "GO:2"]])
    test = _frame([["GO:1"], ["GO:7"], ["GO:8"]])

    vocabulary = select_target_terms(train["terms"].tolist(), min_count=1)
    kept, stats = apply_vocabulary(test, vocabulary, min_terms=1)

    assert vocabulary == ["GO:1", "GO:2"]
    assert list(kept["EntryID"]) == ["P0"]
    assert stats["proteins_dropped"] == 2


def test_apply_vocabulary_does_not_mutate_input():
    """The caller keeps the unfiltered splits for coverage reporting."""
    frame = _frame([["GO:1", "GO:9"]])
    original = frame.copy(deep=True)

    apply_vocabulary(frame, ["GO:1"])

    pd.testing.assert_frame_equal(frame, original)
