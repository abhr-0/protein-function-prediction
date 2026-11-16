"""Unit tests for protein sequence feature extraction modules."""

import numpy as np
import scipy.sparse as sp
import pytest

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
from src.data.uniprot import stream_uniprot_fasta
from src.utils.paths import get_uniprot_path


def test_aac_features():
    """Verify Amino Acid Composition calculation and edge cases."""
    seq = "ACDEFGHIKLMNPQRSTVWY"
    aac = compute_sequence_aac(seq)

    assert aac.shape == (20,)
    # Each standard amino acid appears once in a 20-length sequence -> 1/20 = 0.05
    assert np.allclose(aac, 0.05)
    assert np.isclose(aac.sum(), 1.0)

    # Empty sequence
    empty_aac = compute_sequence_aac("")
    assert np.allclose(empty_aac, 0.0)

    # Sequence with unknown amino acid X
    seq_with_x = "AAAAX"
    aac_x = compute_sequence_aac(seq_with_x)
    assert np.isclose(aac_x[STANDARD_AMINO_ACIDS.index("A")], 4 / 5)

    # Batch extraction
    mat = extract_aac_matrix([seq, "AAAA"])
    assert mat.shape == (2, 20)
    assert len(get_aac_feature_names()) == 20


def test_physicochemical_features():
    """Verify physicochemical properties extraction and robustness to non-standard letters."""
    seq = "MKTIIALSYIFCLVFADYKDDDDK"
    feats = compute_sequence_physicochemical(seq)

    assert feats.shape == (6,)
    length, mw, aroma, pi, gravy, instability = feats

    assert length == 24.0
    assert 2000.0 < mw < 4000.0
    assert 0.0 <= aroma <= 1.0
    assert 2.0 <= pi <= 14.0
    assert np.isfinite(gravy)
    assert np.isfinite(instability)

    # Test ambiguous letter does not raise
    feats_x = compute_sequence_physicochemical("MKTIIALSYIFCLVFADYKDDDDKX")
    assert np.isfinite(feats_x).all()

    # Batch extraction
    mat = extract_physicochemical_matrix([seq, "ACDEF"])
    assert mat.shape == (2, 6)
    assert len(get_physicochemical_feature_names()) == 6


def test_kmer_extractor():
    """Verify K-mer extractor generates normalized sparse matrix."""
    seqs = ["MKTIIALSYIFCLVFADYKDDDDK", "MALLALALALALALALALAL"]
    extractor = KmerExtractor(k=2, normalize_frequencies=True, min_df=1)
    extractor.fit(seqs)

    sparse_mat = extractor.transform(seqs)
    assert sp.issparse(sparse_mat)
    assert sparse_mat.shape[0] == 2

    # Check that row sums are ~1.0 for valid sequences
    row_sums = np.array(sparse_mat.sum(axis=1)).flatten()
    assert np.allclose(row_sums, 1.0)

    names = extractor.get_feature_names_out()
    assert len(names) == sparse_mat.shape[1]
    assert all(name.startswith("kmer_2_") for name in names)


def test_protein_feature_pipeline():
    """Verify unified pipeline fits, transforms, and produces valid finite matrices."""
    train_seqs = [
        "MKTIIALSYIFCLVFADYKDDDDK",
        "MALLALALALALALALALAL",
        "ACDEFGHIKLMNPQRSTVWY",
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
    ]
    test_seqs = ["MKTIIALSYIFCLVFADYKDDDDK", "ACDEFGHIKLMNPQRSTVWY"]

    pipeline = ProteinFeaturePipeline(
        include_aac=True,
        include_physicochemical=True,
        include_kmers=True,
        k=2,
        kmer_min_df=1,
        scale_physicochemical=True,
        return_sparse=True,
    )

    X_train = pipeline.fit_transform(train_seqs)
    assert sp.issparse(X_train)
    assert X_train.shape[0] == len(train_seqs)

    feature_names = pipeline.get_feature_names_out()
    assert len(feature_names) == X_train.shape[1]

    # Transform unseen test set
    X_test = pipeline.transform(test_seqs)
    assert X_test.shape[0] == len(test_seqs)
    assert X_test.shape[1] == X_train.shape[1]

    # Dense option
    pipeline_dense = ProteinFeaturePipeline(
        include_aac=True,
        include_physicochemical=False,
        include_kmers=False,
        return_sparse=False,
    )
    X_dense = pipeline_dense.fit_transform(train_seqs)
    assert isinstance(X_dense, np.ndarray)
    assert X_dense.shape == (len(train_seqs), 20)


def test_pipeline_on_real_uniprot_records():
    """Smoke test pipeline on real sequences streamed from uniprot_sprot.fasta.gz."""
    uniprot_path = get_uniprot_path()
    if not uniprot_path.exists():
        pytest.skip("UniProt file not found.")

    records = list(stream_uniprot_fasta(uniprot_path, max_records=10))
    sequences = [r.sequence for r in records]

    pipeline = ProteinFeaturePipeline(
        include_aac=True,
        include_physicochemical=True,
        include_kmers=True,
        k=2,
        kmer_min_df=1,
    )
    X = pipeline.fit_transform(sequences)

    assert sp.issparse(X)
    assert X.shape[0] == 10
    # Ensure no NaN or Inf in the extracted features
    assert np.all(np.isfinite(X.data))


def test_kmer_vocabulary_is_fitted_on_training_sequences_only():
    """Unseen k-mers must not create feature columns after fitting.

    The k-mer vocabulary and the physicochemical scaler are fitted on training
    sequences; transforming unseen sequences may only reuse those columns.
    """
    pipeline = ProteinFeaturePipeline(
        include_aac=False,
        include_physicochemical=False,
        include_kmers=True,
        k=2,
        kmer_min_df=1,
        return_sparse=True,
    )
    X_train = pipeline.fit_transform(["ACDE", "ACDF"])
    n_features = X_train.shape[1]

    X_unseen = pipeline.transform(["WWWW"])

    assert X_unseen.shape[1] == n_features
    assert X_unseen.nnz == 0, "a sequence of unknown k-mers must be all zeros"
