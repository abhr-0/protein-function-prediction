"""Amino Acid Composition (AAC) feature extraction for protein sequences."""

from typing import Sequence

import numpy as np

# 20 standard IUPAC amino acid single-letter codes
STANDARD_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AMINO_ACID_INDEX = {aa: i for i, aa in enumerate(STANDARD_AMINO_ACIDS)}


def compute_sequence_aac(sequence: str) -> np.ndarray:
    """Compute relative amino acid composition (frequencies) for a single sequence.

    Args:
        sequence: Raw amino acid sequence string.

    Returns:
        1D numpy array of length 20 representing the fraction of each standard amino acid.
        Frequencies sum to 1.0 (or 0.0 if the sequence has no standard amino acids).
    """
    seq_clean = sequence.upper()
    length = len(seq_clean)
    if length == 0:
        return np.zeros(20, dtype=np.float32)

    counts = np.zeros(20, dtype=np.float32)
    for char in seq_clean:
        idx = AMINO_ACID_INDEX.get(char)
        if idx is not None:
            counts[idx] += 1.0

    return counts / length


def extract_aac_matrix(sequences: Sequence[str]) -> np.ndarray:
    """Extract a 2D matrix of AAC features for a collection of sequences.

    Args:
        sequences: Iterable of amino acid sequence strings.

    Returns:
        2D numpy array of shape (N, 20) with dtype float32.
    """
    return np.array([compute_sequence_aac(seq) for seq in sequences], dtype=np.float32)


def get_aac_feature_names() -> list[str]:
    """Return feature names for AAC columns ('aac_A', 'aac_C', etc.)."""
    return [f"aac_{aa}" for aa in STANDARD_AMINO_ACIDS]
