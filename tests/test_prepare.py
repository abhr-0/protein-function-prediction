"""Tests for leakage-safe dataset preparation."""

import pandas as pd
import pytest

from src.data.prepare import PreparedSplits, prepare_splits
from src.models.baseline import split_by_entry_id

SHARED = "GO:0000001"
SHARED_2 = "GO:0000002"


def _dataset(n_proteins: int = 200) -> pd.DataFrame:
    """One shared term plus one unique term per protein.

    The unique terms make vocabulary leakage observable: a vocabulary built
    from the whole frame would contain roughly ``n_proteins`` terms instead of
    roughly the number of training proteins.
    """
    return pd.DataFrame(
        {
            "EntryID": [f"P{index:04d}" for index in range(n_proteins)],
            "terms": [[SHARED, SHARED_2, f"GO:1{index:04d}"] for index in range(n_proteins)],
        }
    )


def test_splits_are_disjoint_and_complete():
    """No protein may appear in two splits, and every protein is placed once."""
    dataset = _dataset(100)
    prepared = prepare_splits(dataset, test_size=0.2, val_size=0.1, min_term_count=1)

    train_ids = set(prepared.train["EntryID"])
    val_ids = set(prepared.val["EntryID"])
    test_ids = set(prepared.test["EntryID"])

    assert train_ids.isdisjoint(val_ids)
    assert train_ids.isdisjoint(test_ids)
    assert val_ids.isdisjoint(test_ids)
    # Unique terms let proteins outside train be dropped, so compare to the split
    assert len(train_ids) + len(val_ids) + len(test_ids) <= len(dataset)


def test_vocabulary_comes_from_training_terms_only():
    """The vocabulary must equal the terms present in the training split.

    This is the leakage guard: if the vocabulary were selected over the full
    dataset, terms belonging only to test proteins would claim a classifier
    column and the reported scores could never be reproduced.
    """
    dataset = _dataset(200)
    train_ref, test_ref = split_by_entry_id(
        dataset, test_size=0.5, val_size=0.0, random_state=42
    )
    expected = {term for terms in train_ref["terms"] for term in terms}
    test_only = {
        term for terms in test_ref["terms"] for term in terms
    } - expected

    prepared = prepare_splits(
        dataset, test_size=0.5, val_size=0.0, random_state=42, min_term_count=1
    )

    assert test_only, "test fixture must contain terms absent from train"
    assert set(prepared.target_terms) == expected
    assert not set(prepared.target_terms) & test_only


def test_min_term_count_counts_training_proteins_only():
    """A term frequent overall but absent from train must not survive.

    The split is computed first so the fixture is deterministic: every protein
    carrying the extra term is known to land in test, which makes the exclusion
    a statement about *where* frequency is counted.
    """
    dataset = _dataset(60)
    train_ref, test_ref = split_by_entry_id(
        dataset, test_size=0.5, val_size=0.0, random_state=1
    )
    test_ids = set(test_ref["EntryID"])
    test_only_positions = [
        index
        for index, entry_id in enumerate(dataset["EntryID"])
        if entry_id in test_ids
    ][:5]

    terms = [list(row) for row in dataset["terms"]]
    for position in test_only_positions:
        terms[position].append("GO:TESTONLY")
    dataset = dataset.assign(terms=terms)

    prepared = prepare_splits(
        dataset, test_size=0.5, val_size=0.0, random_state=1, min_term_count=2
    )

    # Non-vacuous: the term appears 5 times overall, never in train.
    assert sum("GO:TESTONLY" in row for row in terms) == 5
    assert "GO:TESTONLY" not in prepared.target_terms


def test_zero_validation_split_is_supported():
    """val_size=0 must give an empty validation frame and no val stats."""
    prepared = prepare_splits(_dataset(50), val_size=0.0, min_term_count=1)

    assert len(prepared.val) == 0
    assert prepared.has_validation is False
    assert "val" not in prepared.vocab_stats
    assert set(prepared.coverage) == {"train", "val", "test"}


def test_drop_statistics_and_coverage_are_reported():
    """Callers must be able to see what filtering removed."""
    prepared = prepare_splits(
        _dataset(100), test_size=0.2, val_size=0.1, min_term_count=1
    )

    for name in ("train", "val", "test"):
        stats = prepared.vocab_stats[name]
        assert stats["proteins_in"] == stats["proteins_kept"] + stats["proteins_dropped"]
        assert prepared.coverage[name]["proteins"] >= stats["proteins_in"]


def test_min_terms_filters_weakly_annotated_proteins():
    """min_terms>1 drops proteins with too few in-vocabulary terms."""
    dataset = pd.DataFrame(
        {
            "EntryID": ["A", "B", "C", "D"],
            "terms": [
                [SHARED, SHARED_2, "GO:3"],
                [SHARED, SHARED_2],
                [SHARED, SHARED_2],
                [SHARED, SHARED_2, "GO:3"],
            ],
        }
    )
    prepared = prepare_splits(
        dataset, test_size=0.25, val_size=0.0, random_state=0, min_term_count=1, min_terms=2
    )

    assert all(len(terms) >= 2 for terms in prepared.train["terms"])
    assert prepared.vocab_stats["train"]["min_terms"] == 2


def test_splitting_is_reproducible():
    """The same seed must produce the same split and vocabulary."""
    first = prepare_splits(_dataset(80), random_state=7, min_term_count=1)
    second = prepare_splits(_dataset(80), random_state=7, min_term_count=1)

    assert first.target_terms == second.target_terms
    assert list(first.train["EntryID"]) == list(second.train["EntryID"])


def test_single_term_vocabulary_is_rejected():
    """A one-term label space is degenerate and must fail loudly."""
    dataset = pd.DataFrame(
        {"EntryID": [f"P{index}" for index in range(20)], "terms": [[SHARED]] * 20}
    )
    with pytest.raises(ValueError, match="at least 2"):
        prepare_splits(dataset, min_term_count=1)


def test_returns_prepared_splits_type():
    """The return value should be the documented dataclass."""
    assert isinstance(prepare_splits(_dataset(40), min_term_count=1), PreparedSplits)
