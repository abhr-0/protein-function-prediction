"""Leakage-safe dataset preparation: split, vocabulary, zero-positive filtering.

The **order** of operations in this module is the leakage guard:

1. split proteins by ``EntryID``,
2. choose the target vocabulary from the *training* split only,
3. restrict every split to that vocabulary and drop zero-positive proteins.

Choosing the vocabulary before the split lets terms that appear only in test
proteins claim their own classifier column, which inflates scores and cannot be
reproduced on new data.
"""

from dataclasses import dataclass

import pandas as pd

from src.data.targets import apply_vocabulary, select_target_terms, vocabulary_coverage
from src.models.baseline import split_by_entry_id

# NOTE: ``split_by_entry_id`` lives in ``src.models.baseline`` for historical
# reasons. Splitting is a data concern, so it would fit better in ``src/data``;
# moving it is a separate refactor.


@dataclass(frozen=True)
class PreparedSplits:
    """Result of :func:`prepare_splits`, ready for feature extraction."""

    train: pd.DataFrame
    val: pd.DataFrame
    test: pd.DataFrame
    target_terms: list[str]
    coverage: dict[str, dict]
    vocab_stats: dict[str, dict]

    @property
    def has_validation(self) -> bool:
        """Whether a non-empty validation split was produced."""
        return len(self.val) > 0


def prepare_splits(
    dataset: pd.DataFrame,
    *,
    test_size: float = 0.2,
    val_size: float = 0.1,
    random_state: int = 42,
    min_term_count: int = 20,
    top_terms: int | None = None,
    min_terms: int = 1,
    terms_col: str = "terms",
) -> PreparedSplits:
    """Split an aligned dataset and build a train-only target vocabulary.

    Args:
        dataset: Frame with at least ``EntryID`` and a per-protein term list,
            as produced by :func:`src.data.uniprot.build_annotated_dataset`.
        test_size: Fraction of proteins held out for test.
        val_size: Fraction held out for validation; ``0`` gives a 2-way split.
        random_state: Seed for the split, so the same data and seed reproduce it.
        min_term_count: Minimum number of *training* proteins a term must appear
            in to enter the vocabulary.
        top_terms: Optional cap on vocabulary size (``None`` means no cap).
        min_terms: Minimum number of in-vocabulary terms a protein must keep.
        terms_col: Name of the term-list column.

    Returns:
        :class:`PreparedSplits` with filtered frames, the vocabulary, label-space
        coverage per split, and per-split drop statistics.

    Raises:
        ValueError: If fewer than 2 terms survive, or if filtering empties the
            training or test split.
    """
    if val_size > 0:
        train, val, test = split_by_entry_id(
            dataset,
            test_size=test_size,
            val_size=val_size,
            random_state=random_state,
        )
    else:
        train, test = split_by_entry_id(
            dataset,
            test_size=test_size,
            random_state=random_state,
        )
        val = dataset.iloc[0:0].copy()

    # Vocabulary from the raw training terms, before any protein is dropped:
    # re-selecting after filtering would move the label space underneath us.
    target_terms = select_target_terms(
        train[terms_col].tolist(),
        min_count=min_term_count,
        top_k=top_terms,
    )
    if len(target_terms) < 2:
        raise ValueError(
            f"Only {len(target_terms)} target term(s) survived min_count="
            f"{min_term_count}. A multi-label run needs at least 2, so this "
            "configuration is degenerate. Lower min_term_count or use more proteins."
        )

    coverage = {
        name: vocabulary_coverage(frame[terms_col].tolist(), target_terms)
        for name, frame in (("train", train), ("val", val), ("test", test))
    }

    train, train_stats = apply_vocabulary(
        train, target_terms, min_terms=min_terms, terms_col=terms_col
    )
    vocab_stats: dict[str, dict] = {"train": train_stats}
    if len(val):
        val, vocab_stats["val"] = apply_vocabulary(
            val, target_terms, min_terms=min_terms, terms_col=terms_col
        )
    test, vocab_stats["test"] = apply_vocabulary(
        test, target_terms, min_terms=min_terms, terms_col=terms_col
    )

    if len(train) == 0 or len(test) == 0:
        raise ValueError(
            "No proteins remain after vocabulary filtering; loosen "
            "min_term_count, top_terms, or min_terms."
        )

    return PreparedSplits(
        train=train,
        val=val,
        test=test,
        target_terms=target_terms,
        coverage=coverage,
        vocab_stats=vocab_stats,
    )
