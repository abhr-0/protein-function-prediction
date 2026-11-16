"""Tests for basic multi-label metrics and threshold selection."""

import numpy as np
import scipy.sparse as sp

from src.evaluation.metrics import (
    macro_prf,
    micro_prf,
    prf_from_counts,
    select_threshold_micro_f1,
)


def test_prf_from_counts_handles_zero_denominators():
    """Empty denominator cases should return zero, not raise."""
    assert prf_from_counts(0, 0, 0) == {
        "tp": 0,
        "fp": 0,
        "fn": 0,
        "precision": 0.0,
        "recall": 0.0,
        "f1": 0.0,
    }
    assert prf_from_counts(2, 1, 0)["f1"] == 0.8


def test_micro_prf_hand_computed():
    """Micro counts should aggregate over every label cell."""
    y_true = np.array([[1, 0], [0, 1], [0, 0]])
    y_pred = np.array([[1, 0], [1, 1], [0, 0]])

    metrics = micro_prf(y_true, y_pred)

    assert metrics["tp"] == 2
    assert metrics["fp"] == 1
    assert metrics["fn"] == 0
    assert np.isclose(metrics["precision"], 2 / 3)
    assert metrics["recall"] == 1.0
    assert np.isclose(metrics["f1"], 0.8)


def test_macro_prf_skips_empty_terms():
    """Terms with no truth and no predictions should not affect the average."""
    y_true = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]])
    y_pred = np.array([[1, 0, 0], [1, 1, 0], [0, 0, 0]])

    metrics = macro_prf(y_true, y_pred)

    assert metrics["terms"] == 2
    assert np.isclose(metrics["precision"], 0.75)
    assert np.isclose(metrics["recall"], 1.0)
    assert np.isclose(metrics["f1"], (1.0 + 2 / 3) / 2)


def test_metrics_accept_sparse_inputs():
    """Sparse target matrices from the pipeline should work directly."""
    y_true = sp.csr_matrix([[1, 0], [0, 1]])
    y_pred = np.array([[1, 0], [0, 1]])

    metrics = micro_prf(y_true, y_pred)
    assert metrics["f1"] == 1.0


def test_macro_prf_empty_predictions():
    """No predictions should yield zero metrics rather than NaN."""
    metrics = macro_prf(np.zeros((2, 2)), np.zeros((2, 2)))
    assert metrics == {"terms": 0, "precision": 0.0, "recall": 0.0, "f1": 0.0}


def test_select_threshold_prefers_tie_closest_to_half():
    """Tied thresholds should resolve deterministically toward 0.5."""
    y_true = np.array([[1, 0], [0, 1]])
    y_score = np.array([[0.9, 0.1], [0.4, 0.6]])

    threshold, f1 = select_threshold_micro_f1(y_true, y_score)

    assert f1 == 1.0
    assert np.isclose(threshold, 0.5)


def test_select_threshold_without_positives():
    """A label matrix with no positives cannot tune a threshold."""
    threshold, f1 = select_threshold_micro_f1(np.zeros((2, 2)), np.zeros((2, 2)))
    assert threshold == 0.5
    assert f1 == 0.0


def test_select_threshold_rejects_shape_mismatch():
    """Mismatched matrices should fail loudly."""
    try:
        select_threshold_micro_f1(np.ones((2, 2)), np.ones((3, 2)))
    except ValueError as exc:
        assert "Shape mismatch" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("Expected ValueError")
