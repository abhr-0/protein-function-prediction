"""Baseline and advanced prediction models."""

from src.models.baseline import (
    BaselineProteinClassifier,
    build_target_matrix,
    load_baseline_model,
    save_baseline_model,
    split_by_entry_id,
    split_entry_ids,
    target_matrix_to_terms,
)

__all__ = [
    "BaselineProteinClassifier",
    "build_target_matrix",
    "split_by_entry_id",
    "split_entry_ids",
    "target_matrix_to_terms",
    "save_baseline_model",
    "load_baseline_model",
]
