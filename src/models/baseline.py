"""Multi-label baseline model for protein function prediction using Binary Relevance."""

import warnings
from collections.abc import Collection, Sequence
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.base import BaseEstimator, ClassifierMixin, clone
from sklearn.linear_model import LogisticRegression
from sklearn.multiclass import OneVsRestClassifier
from sklearn.utils.validation import check_is_fitted

from src.utils.paths import get_artifacts_dir


def build_target_matrix(
    protein_terms: Sequence[Collection[str]],
    target_terms: Sequence[str],
    return_sparse: bool = True,
    dtype: type = np.int8,
) -> sp.csr_matrix | np.ndarray:
    """Convert per-protein lists of GO terms into a binary indicator matrix.

    Args:
        protein_terms: Sequence of collections of GO term strings for each protein.
        target_terms: Target GO term strings defining the columns.
        return_sparse: If True, returns a scipy.sparse.csr_matrix. Otherwise, dense ndarray.
        dtype: Data type for matrix elements (default np.int8).

    Returns:
        Binary matrix of shape (n_proteins, n_target_terms).
    """
    term_to_idx = {term: idx for idx, term in enumerate(target_terms)}
    rows: list[int] = []
    cols: list[int] = []

    for row_idx, terms in enumerate(protein_terms):
        for term in set(terms):
            col_idx = term_to_idx.get(term)
            if col_idx is not None:
                rows.append(row_idx)
                cols.append(col_idx)

    n_samples = len(protein_terms)
    n_terms = len(target_terms)
    data = np.ones(len(rows), dtype=dtype)
    mat = sp.csr_matrix((data, (rows, cols)), shape=(n_samples, n_terms), dtype=dtype)

    if not return_sparse:
        return mat.toarray()
    return mat


def target_matrix_to_terms(
    matrix: sp.csr_matrix | np.ndarray,
    target_terms: Sequence[str],
    threshold: float = 0.5,
) -> list[list[str]]:
    """Decode a binary or probability matrix back into lists of GO terms per protein.

    Args:
        matrix: Indicator matrix or probability matrix of shape (n_proteins, n_terms).
        target_terms: Target GO term strings corresponding to columns.
        threshold: Decision threshold for inclusion if matrix contains probabilities.

    Returns:
        List of lists containing the selected GO term strings for each protein.
    """
    target_terms_list = list(target_terms)

    if sp.issparse(matrix):
        csr = matrix.tocsr()
        results: list[list[str]] = []
        for i in range(csr.shape[0]):
            row = csr.getrow(i)
            cols = row.indices
            if np.issubdtype(matrix.dtype, np.floating):
                valid_cols = cols[row.data >= threshold]
            else:
                valid_cols = cols[row.data > 0]
            results.append([target_terms_list[c] for c in valid_cols])
        return results
    else:
        results = []
        for row in matrix:
            if np.issubdtype(matrix.dtype, np.floating):
                indices = np.where(row >= threshold)[0]
            else:
                indices = np.where(row > 0)[0]
            results.append([target_terms_list[c] for c in indices])
        return results


def split_entry_ids(
    entry_ids: Sequence[str],
    test_size: float = 0.2,
    val_size: float = 0.0,
    random_state: int = 42,
    shuffle: bool = True,
) -> tuple[list[str], ...]:
    """Split protein entry IDs into train, validation (optional), and test sets.

    Ensures that splitting happens strictly at the protein EntryID level to prevent data leakage.

    Args:
        entry_ids: Sequence of protein accession strings (e.g. UniProt EntryID).
        test_size: Fraction of proteins to allocate to test set.
        val_size: Fraction of proteins to allocate to validation set (0.0 for 2-way split).
        random_state: Seed for reproducibility.
        shuffle: Whether to shuffle before splitting.

    Returns:
        (train_ids, test_ids) if val_size == 0.0, else (train_ids, val_ids, test_ids).
    """
    if not 0.0 < test_size < 1.0:
        raise ValueError(f"test_size must be between 0 and 1, got {test_size}")
    if not 0.0 <= val_size < 1.0:
        raise ValueError(f"val_size must be between 0 and 1, got {val_size}")
    if test_size + val_size >= 1.0:
        raise ValueError(f"Sum of test_size ({test_size}) and val_size ({val_size}) must be < 1.0")

    unique_ids = list(dict.fromkeys(entry_ids))
    n_total = len(unique_ids)
    if n_total == 0:
        raise ValueError("Cannot split an empty sequence of entry IDs.")

    indices = np.arange(n_total)
    if shuffle:
        rng = np.random.default_rng(random_state)
        rng.shuffle(indices)

    n_test = max(1, int(np.round(n_total * test_size))) if test_size > 0 else 0
    n_val = int(np.round(n_total * val_size)) if val_size > 0 else 0
    if n_test + n_val >= n_total:
        n_test = max(1, n_total - n_val - 1)

    n_train = n_total - n_test - n_val

    train_idx = indices[:n_train]
    val_idx = indices[n_train : n_train + n_val]
    test_idx = indices[n_train + n_val :]

    train_ids = [unique_ids[i] for i in train_idx]
    test_ids = [unique_ids[i] for i in test_idx]

    if val_size > 0:
        val_ids = [unique_ids[i] for i in val_idx]
        return train_ids, val_ids, test_ids
    return train_ids, test_ids


def split_by_entry_id(
    data: pd.DataFrame,
    test_size: float = 0.2,
    val_size: float = 0.0,
    random_state: int = 42,
    id_col: str = "EntryID",
    shuffle: bool = True,
) -> tuple[pd.DataFrame, ...]:
    """Split a DataFrame of proteins by unique EntryID to prevent data leakage.

    Args:
        data: DataFrame containing an EntryID column.
        test_size: Proportion of unique proteins in the test set.
        val_size: Proportion of unique proteins in the validation set.
        random_state: Random seed for shuffling.
        id_col: Column name containing protein entry identifiers.
        shuffle: Whether to shuffle proteins before splitting.

    Returns:
        (train_df, test_df) if val_size == 0.0, else (train_df, val_df, test_df).
    """
    if id_col not in data.columns:
        raise KeyError(f"Column '{id_col}' not found in DataFrame.")

    entry_ids = data[id_col].tolist()
    split_results = split_entry_ids(
        entry_ids=entry_ids,
        test_size=test_size,
        val_size=val_size,
        random_state=random_state,
        shuffle=shuffle,
    )

    if val_size > 0:
        train_ids, val_ids, test_ids = split_results
        train_df = data[data[id_col].isin(set(train_ids))].copy().reset_index(drop=True)
        val_df = data[data[id_col].isin(set(val_ids))].copy().reset_index(drop=True)
        test_df = data[data[id_col].isin(set(test_ids))].copy().reset_index(drop=True)
        return train_df, val_df, test_df
    else:
        train_ids, test_ids = split_results
        train_df = data[data[id_col].isin(set(train_ids))].copy().reset_index(drop=True)
        test_df = data[data[id_col].isin(set(test_ids))].copy().reset_index(drop=True)
        return train_df, test_df


class BaselineProteinClassifier(BaseEstimator, ClassifierMixin):
    """Multi-label classifier using Binary Relevance (One-vs-Rest) with calibrated probabilities.

    Wraps a base binary classifier (by default, LogisticRegression with liblinear solver)
    inside OneVsRestClassifier with multi-core parallelization.
    """

    def __init__(
        self,
        base_estimator: BaseEstimator | None = None,
        target_terms: Sequence[str] | None = None,
        n_jobs: int = -1,
        random_state: int = 42,
    ) -> None:
        self.base_estimator = base_estimator
        self.target_terms = target_terms
        self.n_jobs = n_jobs
        self.random_state = random_state

        self.classifier_: OneVsRestClassifier | None = None
        self.target_terms_: list[str] | None = None
        self.n_classes_: int = 0

    def fit(
        self,
        X: sp.csr_matrix | np.ndarray,
        y: sp.csr_matrix | np.ndarray,
        target_terms: Sequence[str] | None = None,
    ) -> "BaselineProteinClassifier":
        """Fit binary relevance classifiers for each target GO term.

        Args:
            X: Protein feature matrix of shape (n_proteins, n_features).
            y: Binary target matrix of shape (n_proteins, n_terms).
            target_terms: Optional list of target GO terms matching the columns of y.

        Returns:
            self
        """
        if X.shape[0] != y.shape[0]:
            raise ValueError(
                f"Feature matrix X has {X.shape[0]} samples but target matrix y has {y.shape[0]}."
            )

        n_classes = y.shape[1] if y.ndim > 1 else 1
        if y.ndim > 1 and n_classes < 2:
            raise ValueError(
                "Multi-label training needs at least 2 target columns, got 1. With a "
                "single column, sklearn treats the task as multiclass and fails on "
                "sparse input (unhashable type: 'csr_matrix'). A one-term vocabulary "
                "is degenerate: lower the term frequency threshold or use more proteins."
            )
        self.n_classes_ = n_classes

        if target_terms is not None:
            self.target_terms_ = list(target_terms)
        elif self.target_terms is not None:
            self.target_terms_ = list(self.target_terms)
        else:
            self.target_terms_ = [f"term_{i}" for i in range(n_classes)]

        if len(self.target_terms_) != self.n_classes_:
            raise ValueError(
                f"Target terms length ({len(self.target_terms_)}) does not match y columns ({self.n_classes_})."
            )

        if self.base_estimator is None:
            estimator = LogisticRegression(
                solver="liblinear",
                C=1.0,
                max_iter=200,
                random_state=self.random_state,
            )
        else:
            estimator = clone(self.base_estimator)

        self.classifier_ = OneVsRestClassifier(estimator=estimator, n_jobs=self.n_jobs)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=r".*present in all training examples.*")
            self.classifier_.fit(X, y)
        return self

    def predict_proba(self, X: sp.csr_matrix | np.ndarray) -> np.ndarray:
        """Predict probability scores for each target GO term.

        Args:
            X: Feature matrix of shape (n_proteins, n_features).

        Returns:
            Probability matrix of shape (n_proteins, n_terms) with values in [0.0, 1.0].
        """
        check_is_fitted(self, "classifier_")
        assert self.classifier_ is not None

        probs = self.classifier_.predict_proba(X)

        # Edge case: if single class was passed, OneVsRestClassifier outputs (N, 2)
        if self.n_classes_ == 1 and probs.ndim == 2 and probs.shape[1] == 2:
            probs = probs[:, [1]]

        return probs

    def predict(
        self,
        X: sp.csr_matrix | np.ndarray,
        threshold: float = 0.5,
    ) -> np.ndarray:
        """Predict binary indicator matrix at a given probability threshold.

        Args:
            X: Feature matrix of shape (n_proteins, n_features).
            threshold: Probability decision threshold (default 0.5).

        Returns:
            Binary indicator matrix of shape (n_proteins, n_terms).
        """
        probs = self.predict_proba(X)
        return (probs >= threshold).astype(np.int8)

    def predict_top_k(
        self,
        X: sp.csr_matrix | np.ndarray,
        k: int = 10,
        threshold: float = 0.0,
    ) -> list[list[tuple[str, float]]]:
        """Predict top-k most confident GO terms for each protein.

        Args:
            X: Feature matrix of shape (n_proteins, n_features).
            k: Maximum number of terms to return per protein.
            threshold: Minimum probability required for inclusion.

        Returns:
            List of lists of (term_id, probability) sorted in descending order of confidence.
        """
        probs = self.predict_proba(X)
        terms = self.target_terms_ or [f"term_{i}" for i in range(probs.shape[1])]

        results: list[list[tuple[str, float]]] = []
        for row in probs:
            sorted_indices = np.argsort(row)[::-1]
            top_terms = [
                (terms[idx], float(row[idx]))
                for idx in sorted_indices
                if row[idx] >= threshold
            ][:k]
            results.append(top_terms)

        return results

    def save(self, filepath: str | Path | None = None) -> Path:
        """Serialize and save the trained model to disk.

        Args:
            filepath: Destination file path. If None, saves to data/artifacts/baseline_model.joblib.

        Returns:
            Path to the saved artifact.
        """
        check_is_fitted(self, "classifier_")
        if filepath is None:
            dest = get_artifacts_dir() / "baseline_model.joblib"
        else:
            dest = Path(filepath)

        dest.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(self, dest)
        return dest

    @classmethod
    def load(cls, filepath: str | Path) -> "BaselineProteinClassifier":
        """Load a saved model from disk.

        Args:
            filepath: Path to the serialized .joblib file.

        Returns:
            Loaded BaselineProteinClassifier instance.
        """
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"Model file not found: {path}")

        loaded = joblib.load(path)
        if not isinstance(loaded, cls):
            raise TypeError(f"Loaded object is {type(loaded)}, expected {cls}.")
        return loaded


def save_baseline_model(
    model: BaselineProteinClassifier,
    filepath: str | Path | None = None,
) -> Path:
    """Save a BaselineProteinClassifier instance."""
    return model.save(filepath)


def load_baseline_model(filepath: str | Path) -> BaselineProteinClassifier:
    """Load a BaselineProteinClassifier instance."""
    return BaselineProteinClassifier.load(filepath)
