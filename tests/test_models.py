"""Tests for multi-label baseline models, target matrix builder, and split utilities."""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from src.features.pipeline import ProteinFeaturePipeline
from src.models.baseline import (
    BaselineProteinClassifier,
    build_target_matrix,
    load_baseline_model,
    save_baseline_model,
    split_by_entry_id,
    split_entry_ids,
    target_matrix_to_terms,
)


def test_build_target_matrix_sparse():
    protein_terms = [
        ["GO:0005515", "GO:0005524"],
        ["GO:0005524"],
        ["GO:9999999"],  # Not in target_terms
        [],
    ]
    target_terms = ["GO:0005515", "GO:0005524", "GO:0003677"]

    y = build_target_matrix(protein_terms, target_terms, return_sparse=True)
    assert sp.issparse(y)
    assert y.shape == (4, 3)

    dense_y = y.toarray()
    assert (dense_y[0] == [1, 1, 0]).all()
    assert (dense_y[1] == [0, 1, 0]).all()
    assert (dense_y[2] == [0, 0, 0]).all()
    assert (dense_y[3] == [0, 0, 0]).all()


def test_build_target_matrix_dense():
    protein_terms = [["GO:0001", "GO:0002"], ["GO:0002"]]
    target_terms = ["GO:0001", "GO:0002"]

    y = build_target_matrix(protein_terms, target_terms, return_sparse=False)
    assert isinstance(y, np.ndarray)
    assert y.shape == (2, 2)
    assert (y[0] == [1, 1]).all()
    assert (y[1] == [0, 1]).all()


def test_target_matrix_to_terms():
    target_terms = ["GO:0001", "GO:0002", "GO:0003"]

    # Test round trip with sparse binary matrix
    y_sparse = sp.csr_matrix([[1, 0, 1], [0, 1, 0], [0, 0, 0]])
    decoded = target_matrix_to_terms(y_sparse, target_terms)
    assert decoded == [["GO:0001", "GO:0003"], ["GO:0002"], []]

    # Test with dense probability matrix and threshold
    y_prob = np.array([[0.8, 0.2, 0.6], [0.1, 0.9, 0.4]])
    decoded_prob = target_matrix_to_terms(y_prob, target_terms, threshold=0.5)
    assert decoded_prob == [["GO:0001", "GO:0003"], ["GO:0002"]]


def test_split_entry_ids():
    entry_ids = [f"P{i:03d}" for i in range(100)]

    # 2-way split
    train_ids, test_ids = split_entry_ids(entry_ids, test_size=0.2, random_state=42)
    assert len(train_ids) == 80
    assert len(test_ids) == 20
    assert set(train_ids).isdisjoint(set(test_ids))

    # 3-way split
    train_ids, val_ids, test_ids = split_entry_ids(
        entry_ids, test_size=0.2, val_size=0.1, random_state=42
    )
    assert len(train_ids) == 70
    assert len(val_ids) == 10
    assert len(test_ids) == 20
    assert set(train_ids).isdisjoint(set(val_ids))
    assert set(train_ids).isdisjoint(set(test_ids))
    assert set(val_ids).isdisjoint(set(test_ids))

    # Error handling for invalid splits
    with pytest.raises(ValueError, match="Sum of test_size"):
        split_entry_ids(entry_ids, test_size=0.6, val_size=0.5)


def test_split_by_entry_id():
    df = pd.DataFrame(
        {
            "EntryID": [f"P{i}" for i in range(20)],
            "sequence": ["MKTIIALSYIFCLVFA"] * 20,
            "terms": [["GO:0001"]] * 20,
        }
    )

    train_df, test_df = split_by_entry_id(df, test_size=0.25, random_state=42)
    assert len(train_df) == 15
    assert len(test_df) == 5
    assert set(train_df["EntryID"]).isdisjoint(set(test_df["EntryID"]))

    # Test missing column
    with pytest.raises(KeyError, match="Column 'non_existent' not found"):
        split_by_entry_id(df, id_col="non_existent")


def test_baseline_protein_classifier_fit_predict():
    np.random.seed(42)
    n_samples = 40
    n_features = 20
    n_classes = 5

    X = sp.csr_matrix(np.random.randn(n_samples, n_features))
    y = sp.csr_matrix(np.random.randint(0, 2, size=(n_samples, n_classes), dtype=np.int8))
    target_terms = [f"GO:{i:07d}" for i in range(n_classes)]

    clf = BaselineProteinClassifier(target_terms=target_terms, n_jobs=1, random_state=42)
    clf.fit(X, y)

    # Predict probabilities
    probs = clf.predict_proba(X)
    assert probs.shape == (n_samples, n_classes)
    assert (probs >= 0.0).all() and (probs <= 1.0).all()

    # Predict binary indicators
    preds = clf.predict(X, threshold=0.5)
    assert preds.shape == (n_samples, n_classes)
    assert set(np.unique(preds)).issubset({0, 1})


def test_baseline_protein_classifier_predict_top_k():
    n_samples = 5
    n_classes = 4
    X = sp.csr_matrix(np.random.randn(n_samples, 10))
    y = sp.csr_matrix(np.random.randint(0, 2, size=(n_samples, n_classes), dtype=np.int8))
    target_terms = ["GO:0001", "GO:0002", "GO:0003", "GO:0004"]

    clf = BaselineProteinClassifier(target_terms=target_terms, n_jobs=1)
    clf.fit(X, y)

    top_k = clf.predict_top_k(X, k=2, threshold=0.0)
    assert len(top_k) == n_samples
    for row in top_k:
        assert len(row) <= 2
        # Verify probabilities are sorted descending
        if len(row) == 2:
            assert row[0][1] >= row[1][1]
        assert row[0][0] in target_terms


def test_baseline_protein_classifier_save_and_load():
    X = sp.csr_matrix(np.random.randn(20, 10))
    y = sp.csr_matrix(np.random.randint(0, 2, size=(20, 3), dtype=np.int8))
    target_terms = ["GO:0001", "GO:0002", "GO:0003"]

    clf = BaselineProteinClassifier(target_terms=target_terms, n_jobs=1)
    clf.fit(X, y)
    probs_orig = clf.predict_proba(X)

    with tempfile.TemporaryDirectory() as tmp_dir:
        save_path = Path(tmp_dir) / "test_model.joblib"
        saved_path = save_baseline_model(clf, save_path)
        assert saved_path.exists()

        loaded_clf = load_baseline_model(saved_path)
        probs_loaded = loaded_clf.predict_proba(X)

        assert np.allclose(probs_orig, probs_loaded)
        assert loaded_clf.target_terms_ == target_terms


def test_baseline_pipeline_integration():
    """End-to-end integration test: sequences -> features -> targets -> baseline fit & predict."""
    sequences = [
        "MKTIIALSYIFCLVFA",
        "ACDEFGHIKLMNPQRSTVWY",
        "MKTIIALSYIFCLVFADYKDDDDK",
        "KKKKKKKKKKKKKKKK",
        "WFWFWFWFWFWFWFWF",
        "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY",
        "MKTIIALSYIFCLVFAWFWF",
        "ACDEFACDEFACDEF",
    ]
    protein_terms = [
        ["GO:0005515", "GO:0005524"],
        ["GO:0005524"],
        ["GO:0005515"],
        ["GO:0003677"],
        ["GO:0003677", "GO:0005515"],
        ["GO:0005524"],
        ["GO:0005515"],
        ["GO:0005524", "GO:0003677"],
    ]
    target_terms = ["GO:0005515", "GO:0005524", "GO:0003677"]

    # 1. Feature extraction
    feature_pipeline = ProteinFeaturePipeline(k=2, kmer_min_df=1, return_sparse=True)
    X = feature_pipeline.fit_transform(sequences)
    assert sp.issparse(X)
    assert X.shape[0] == len(sequences)

    # 2. Target matrix
    y = build_target_matrix(protein_terms, target_terms, return_sparse=True)
    assert y.shape == (len(sequences), len(target_terms))

    # 3. Model training
    model = BaselineProteinClassifier(target_terms=target_terms, n_jobs=1)
    model.fit(X, y)

    # 4. Predict on a new sequence
    new_seq = ["MKTIIALSYIFCLVFAACDEF"]
    X_new = feature_pipeline.transform(new_seq)
    probs_new = model.predict_proba(X_new)
    assert probs_new.shape == (1, len(target_terms))
    assert (probs_new >= 0.0).all() and (probs_new <= 1.0).all()

    top_predictions = model.predict_top_k(X_new, k=2)
    assert len(top_predictions) == 1
    assert len(top_predictions[0]) == 2
