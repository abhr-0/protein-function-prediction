"""Unified protein sequence feature extraction pipeline."""

from typing import Sequence

import numpy as np
import scipy.sparse as sp
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler

from src.features.aac import extract_aac_matrix, get_aac_feature_names
from src.features.kmers import KmerExtractor
from src.features.physicochemical import (
    extract_physicochemical_matrix,
    get_physicochemical_feature_names,
)


class ProteinFeaturePipeline(BaseEstimator, TransformerMixin):
    """Combines AAC, physicochemical descriptors, and k-mer frequencies into a unified matrix."""

    def __init__(
        self,
        include_aac: bool = True,
        include_physicochemical: bool = True,
        include_kmers: bool = True,
        k: int = 2,
        kmer_min_df: int | float = 2,
        scale_physicochemical: bool = True,
        return_sparse: bool = True,
    ) -> None:
        self.include_aac = include_aac
        self.include_physicochemical = include_physicochemical
        self.include_kmers = include_kmers
        self.k = k
        self.kmer_min_df = kmer_min_df
        self.scale_physicochemical = scale_physicochemical
        self.return_sparse = return_sparse

        self._kmer_extractor: KmerExtractor | None = None
        self._scaler: StandardScaler | None = None
        self._feature_names: list[str] = []

    def fit(self, sequences: Sequence[str], y: None = None) -> "ProteinFeaturePipeline":
        """Fit vocabulary for k-mers and scalers for physicochemical features."""
        feature_names = []

        if self.include_aac:
            feature_names.extend(get_aac_feature_names())

        if self.include_physicochemical:
            phys_names = get_physicochemical_feature_names()
            feature_names.extend(phys_names)
            if self.scale_physicochemical:
                phys_matrix = extract_physicochemical_matrix(sequences)
                self._scaler = StandardScaler()
                self._scaler.fit(phys_matrix)

        if self.include_kmers:
            self._kmer_extractor = KmerExtractor(
                k=self.k,
                normalize_frequencies=True,
                min_df=self.kmer_min_df,
            )
            self._kmer_extractor.fit(sequences)
            feature_names.extend(self._kmer_extractor.get_feature_names_out())

        self._feature_names = feature_names
        return self

    def transform(self, sequences: Sequence[str]) -> sp.csr_matrix | np.ndarray:
        """Transform a collection of protein sequences into a feature matrix."""
        matrices_to_stack = []

        if self.include_aac:
            aac_mat = extract_aac_matrix(sequences)
            matrices_to_stack.append(sp.csr_matrix(aac_mat))

        if self.include_physicochemical:
            phys_mat = extract_physicochemical_matrix(sequences)
            if self.scale_physicochemical and self._scaler is not None:
                phys_mat = self._scaler.transform(phys_mat)
            matrices_to_stack.append(sp.csr_matrix(phys_mat))

        if self.include_kmers and self._kmer_extractor is not None:
            kmer_mat = self._kmer_extractor.transform(sequences)
            matrices_to_stack.append(kmer_mat)

        if not matrices_to_stack:
            raise ValueError("At least one feature family (aac, physicochemical, kmers) must be enabled.")

        combined = sp.hstack(matrices_to_stack).tocsr()

        if not self.return_sparse:
            return combined.toarray()
        return combined

    def get_feature_names_out(self) -> list[str]:
        """Return the combined list of feature column names."""
        return list(self._feature_names)
