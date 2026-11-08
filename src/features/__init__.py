"""Protein sequence feature extraction modules."""

from src.features.aac import (
    STANDARD_AMINO_ACIDS,
    compute_sequence_aac,
    extract_aac_matrix,
    get_aac_feature_names,
)
from src.features.kmers import KmerExtractor
from src.features.physicochemical import (
    compute_sequence_physicochemical,
    extract_physicochemical_matrix,
    get_physicochemical_feature_names,
)
from src.features.pipeline import ProteinFeaturePipeline

__all__ = [
    "STANDARD_AMINO_ACIDS",
    "compute_sequence_aac",
    "extract_aac_matrix",
    "get_aac_feature_names",
    "compute_sequence_physicochemical",
    "extract_physicochemical_matrix",
    "get_physicochemical_feature_names",
    "KmerExtractor",
    "ProteinFeaturePipeline",
]
