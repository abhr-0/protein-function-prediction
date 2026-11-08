"""K-mer feature extraction for protein sequences using sparse matrices."""

from typing import Sequence

import scipy.sparse as sp
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize


class KmerExtractor(BaseEstimator, TransformerMixin):
    """Extracts k-mer counts or normalized frequencies from amino acid sequences.

    Produces memory-efficient `scipy.sparse.csr_matrix` outputs.
    """

    def __init__(
        self,
        k: int = 2,
        normalize_frequencies: bool = True,
        min_df: int | float = 1,
        max_features: int | None = None,
    ) -> None:
        self.k = k
        self.normalize_frequencies = normalize_frequencies
        self.min_df = min_df
        self.max_features = max_features
        self._vectorizer: CountVectorizer | None = None

    def fit(self, sequences: Sequence[str], y: None = None) -> "KmerExtractor":
        """Fit vocabulary from training sequences."""
        self._vectorizer = CountVectorizer(
            analyzer="char",
            ngram_range=(self.k, self.k),
            lowercase=False,
            min_df=self.min_df,
            max_features=self.max_features,
        )
        # Ensure all sequences are uppercase strings
        cleaned = [str(s).upper() for s in sequences]
        self._vectorizer.fit(cleaned)
        return self

    def transform(self, sequences: Sequence[str]) -> sp.csr_matrix:
        """Transform sequences into sparse k-mer feature matrix."""
        if self._vectorizer is None:
            raise RuntimeError("KmerExtractor has not been fitted yet.")

        cleaned = [str(s).upper() for s in sequences]
        counts = self._vectorizer.transform(cleaned)

        if self.normalize_frequencies:
            return normalize(counts, norm="l1", axis=1)
        return counts

    def get_feature_names_out(self) -> list[str]:
        """Return feature names for all vocabulary k-mers (e.g. 'kmer_2_AL')."""
        if self._vectorizer is None:
            raise RuntimeError("KmerExtractor has not been fitted yet.")
        raw_names = self._vectorizer.get_feature_names_out()
        return [f"kmer_{self.k}_{name}" for name in raw_names]
