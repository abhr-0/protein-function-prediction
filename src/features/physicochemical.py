"""Physicochemical feature extraction using Biopython ProtParam."""

from typing import Sequence

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from src.features.aac import STANDARD_AMINO_ACIDS

PHYSICOCHEMICAL_FEATURE_NAMES = [
    "length",
    "molecular_weight",
    "aromaticity",
    "isoelectric_point",
    "gravy",
    "instability_index",
]


def compute_sequence_physicochemical(sequence: str) -> np.ndarray:
    """Compute physicochemical descriptors for a single protein sequence.

    Strips ambiguous/non-standard characters (X, B, Z, U) before ProtParam analysis.

    Returns:
        1D numpy array of 6 floats: [length, molecular_weight, aromaticity, pI, gravy, instability].
    """
    raw_len = len(sequence)
    if raw_len == 0:
        return np.zeros(6, dtype=np.float32)

    # Clean to unambiguous 20 standard amino acids for Biopython
    clean_seq = "".join(c for c in sequence.upper() if c in STANDARD_AMINO_ACIDS)

    if not clean_seq:
        return np.array([float(raw_len), 0.0, 0.0, 7.0, 0.0, 0.0], dtype=np.float32)

    try:
        analysis = ProteinAnalysis(clean_seq)
        mw = analysis.molecular_weight()
        aroma = analysis.aromaticity()
        pi = analysis.isoelectric_point()
        gravy = analysis.gravy()
        instability = analysis.instability_index()
    except Exception:
        # Fallback in case of edge cases in calculation
        mw = float(raw_len * 110.0)  # average amino acid MW ~110 Da
        aroma = 0.0
        pi = 7.0
        gravy = 0.0
        instability = 0.0

    return np.array([float(raw_len), mw, aroma, pi, gravy, instability], dtype=np.float32)


def extract_physicochemical_matrix(sequences: Sequence[str]) -> np.ndarray:
    """Extract a 2D matrix of physicochemical descriptors for a list of sequences.

    Returns:
        2D numpy array of shape (N, 6) with dtype float32.
    """
    return np.array(
        [compute_sequence_physicochemical(seq) for seq in sequences],
        dtype=np.float32,
    )


def get_physicochemical_feature_names() -> list[str]:
    """Return feature names for physicochemical columns."""
    return list(PHYSICOCHEMICAL_FEATURE_NAMES)
