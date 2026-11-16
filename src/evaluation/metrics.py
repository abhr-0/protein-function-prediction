"""Basic multi-label metrics and threshold selection.

This module starts with threshold-independent confusion-matrix helpers used by
the Phase 1 sanity check. Hierarchy-aware CAFA metrics are added in Phase 2.
"""

import numpy as np
import scipy.sparse as sp


def _to_dense(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def prf_from_counts(tp: int, fp: int, fn: int) -> dict[str, float | int]:
    """Compute precision, recall, and F1 from confusion counts."""
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {
        "tp": int(tp),
        "fp": int(fp),
        "fn": int(fn),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
    }


def micro_prf(y_true, y_pred) -> dict[str, float | int]:
    """Compute micro-averaged precision/recall/F1 across all label cells."""
    true = _to_dense(y_true).astype(bool)
    pred = _to_dense(y_pred).astype(bool)
    if true.shape != pred.shape:
        raise ValueError(f"Shape mismatch: y_true {true.shape} vs y_pred {pred.shape}")

    tp = int(np.sum(true & pred))
    fp = int(np.sum(~true & pred))
    fn = int(np.sum(true & ~pred))
    return prf_from_counts(tp, fp, fn)


def macro_prf(y_true, y_pred) -> dict[str, float | int]:
    """Average precision/recall/F1 over terms with support or predictions.

    Terms with neither ground-truth positives nor predictions are skipped, so a
    truncated vocabulary does not inflate the macro average with empty columns.
    """
    true = _to_dense(y_true).astype(bool)
    pred = _to_dense(y_pred).astype(bool)
    if true.shape != pred.shape:
        raise ValueError(f"Shape mismatch: y_true {true.shape} vs y_pred {pred.shape}")

    tp = np.sum(true & pred, axis=0)
    fp = np.sum(~true & pred, axis=0)
    fn = np.sum(true & ~pred, axis=0)

    valid = (tp + fp + fn) > 0
    if not np.any(valid):
        return {"terms": 0, "precision": 0.0, "recall": 0.0, "f1": 0.0}

    tp = tp[valid].astype(float)
    fp = fp[valid].astype(float)
    fn = fn[valid].astype(float)

    precision = np.divide(tp, tp + fp, out=np.zeros_like(tp), where=(tp + fp) > 0)
    recall = np.divide(tp, tp + fn, out=np.zeros_like(tp), where=(tp + fn) > 0)
    f1 = np.divide(
        2 * precision * recall,
        precision + recall,
        out=np.zeros_like(tp),
        where=(precision + recall) > 0,
    )
    return {
        "terms": int(np.sum(valid)),
        "precision": float(precision.mean()),
        "recall": float(recall.mean()),
        "f1": float(f1.mean()),
    }


def select_threshold_micro_f1(
    y_true,
    y_score,
    thresholds: np.ndarray | None = None,
) -> tuple[float, float]:
    """Choose the probability threshold that maximizes micro-F1.

    Ties are broken by the threshold closest to 0.5, then by the lower
    threshold, so the result is deterministic.
    """
    true = _to_dense(y_true).astype(bool)
    score = _to_dense(y_score)
    if true.shape != score.shape:
        raise ValueError(f"Shape mismatch: y_true {true.shape} vs y_score {score.shape}")

    if thresholds is None:
        thresholds = np.linspace(0.0, 1.0, 101)

    if true.sum() == 0:
        return 0.5, 0.0

    best_f1 = -1.0
    best_thresholds: list[float] = []
    for threshold in thresholds:
        pred = score >= threshold
        f1 = micro_prf(true, pred)["f1"]
        if f1 > best_f1 + 1e-12:
            best_f1 = float(f1)
            best_thresholds = [float(threshold)]
        elif abs(f1 - best_f1) <= 1e-12:
            best_thresholds.append(float(threshold))

    best_threshold = min(best_thresholds, key=lambda value: (abs(value - 0.5), value))
    return best_threshold, best_f1
