#!/usr/bin/env python3
"""Run baseline protein function prediction pipeline on a small subset of proteins."""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

# Ensure project root is in pythonpath
project_root = Path(__file__).resolve().parents[1]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.data.uniprot import create_annotated_dataset
from src.features.pipeline import ProteinFeaturePipeline
from src.models.baseline import (
    BaselineProteinClassifier,
    build_target_matrix,
    split_by_entry_id,
    target_matrix_to_terms,
)
from src.utils.paths import (
    get_artifacts_dir,
    get_train_terms_path,
    get_uniprot_path,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train and evaluate baseline classifier on a subset of UniProt Swiss-Prot proteins."
    )
    parser.add_argument(
        "--n-proteins",
        type=int,
        default=200,
        help="Number of annotated proteins to extract (0 or negative for full dataset, default: 200)",
    )
    parser.add_argument(
        "--top-terms",
        type=int,
        default=50,
        help="Number of most frequent target GO terms to predict (default: 50)",
    )
    parser.add_argument(
        "--aspect",
        type=str,
        default="F",
        help="Gene Ontology aspect: 'F' (Molecular Function), 'P' (Biological Process), 'C' (Cellular Component)",
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.2,
        help="Fraction of proteins for test set (default: 0.2)",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed for data splitting and model training (default: 42)",
    )
    args = parser.parse_args()

    max_proteins = args.n_proteins if args.n_proteins > 0 else None
    dataset_label = f"{args.n_proteins} proteins" if max_proteins else "FULL DATASET (all Swiss-Prot MF proteins)"

    print("=" * 70)
    print("PROTEIN FUNCTION PREDICTION — BASELINE RUN")
    print("=" * 70)
    print(f"Dataset Target: {dataset_label}")
    print(f"Target Terms:   Top {args.top_terms} GO terms (Aspect: {args.aspect})")
    print(f"Test Split:     {args.test_size * 100:.0f}%")
    print("-" * 70)

    # 1. Load Data
    t0 = time.time()
    print("1. Ingesting sequences & annotations from disk...")
    df, target_terms = create_annotated_dataset(
        fasta_path=get_uniprot_path(),
        terms_path=get_train_terms_path(),
        aspect=args.aspect,
        max_proteins=max_proteins,
        top_k_terms=args.top_terms,
    )
    t_load = time.time() - t0
    print(
        f"   ✓ Loaded {len(df)} matched proteins with {len(target_terms)} target GO terms in {t_load:.2f}s."
    )

    # 2. Split by EntryID to prevent leakage
    print("2. Splitting dataset by EntryID...")
    train_df, test_df = split_by_entry_id(
        df,
        test_size=args.test_size,
        random_state=args.random_state,
    )
    print(f"   ✓ Train set: {len(train_df)} proteins")
    print(f"   ✓ Test set:  {len(test_df)} proteins")

    # 3. Extract Features
    t0 = time.time()
    print("3. Extracting features (AAC + Physicochemical + 2-mers)...")
    feature_pipeline = ProteinFeaturePipeline(
        include_aac=True,
        include_physicochemical=True,
        include_kmers=True,
        k=2,
        kmer_min_df=2,
        scale_physicochemical=True,
        return_sparse=True,
    )
    X_train = feature_pipeline.fit_transform(train_df["sequence"].tolist())
    X_test = feature_pipeline.transform(test_df["sequence"].tolist())
    t_feat = time.time() - t0
    print(
        f"   ✓ Features extracted in {t_feat:.2f}s | Shape: Train {X_train.shape}, Test {X_test.shape}"
    )
    print(f"   ✓ Total sequence descriptors: {X_train.shape[1]}")

    # 4. Build Multi-Label Targets
    print("4. Building binary target matrices...")
    y_train = build_target_matrix(train_df["terms"].tolist(), target_terms, return_sparse=True)
    y_test = build_target_matrix(test_df["terms"].tolist(), target_terms, return_sparse=True)
    print(f"   ✓ Target matrix shape: {y_train.shape} (sparsity: {100.0 * (1.0 - y_train.nnz / (y_train.shape[0] * y_train.shape[1])):.1f}%)")

    # 5. Train Model
    t0 = time.time()
    print("5. Training BaselineProteinClassifier (One-vs-Rest Logistic Regression)...")
    model = BaselineProteinClassifier(
        target_terms=target_terms,
        n_jobs=-1,
        random_state=args.random_state,
    )
    model.fit(X_train, y_train)
    t_train = time.time() - t0
    print(f"   ✓ Model training completed in {t_train:.2f}s across all CPU cores.")

    # 6. Predict on Test Set
    t0 = time.time()
    print("6. Predicting on test set...")
    test_probs = model.predict_proba(X_test)
    test_preds_binary = model.predict(X_test, threshold=0.5)
    t_pred = time.time() - t0
    print(f"   ✓ Predictions generated in {t_pred:.3f}s.")

    # 7. Model Persistence
    tag = f"{args.n_proteins}" if max_proteins else "full"
    artifact_path = get_artifacts_dir() / f"baseline_model_{tag}.joblib"
    saved_path = model.save(artifact_path)
    print(f"7. Serialized model saved to: {saved_path}")

    # 8. Inspection / Sample Predictions
    print("-" * 70)
    print("SAMPLE PREDICTIONS (First 3 Test Proteins):")
    print("-" * 70)
    top_k_preds = model.predict_top_k(X_test[:3], k=5, threshold=0.1)

    for i in range(min(3, len(test_df))):
        entry_id = test_df.iloc[i]["EntryID"]
        name = test_df.iloc[i]["entry_name"]
        org = test_df.iloc[i]["organism"]
        true_terms = test_df.iloc[i]["terms"]
        length = test_df.iloc[i]["length"]

        print(f"\nProtein #{i+1}: {entry_id} ({name}) | Length: {length} aa | Organism: {org}")
        print(f"  • Ground Truth MF Terms ({len(true_terms)}): {', '.join(true_terms[:6])}")
        print("  • Top Predicted Terms (Confidence):")
        if top_k_preds[i]:
            for term, prob in top_k_preds[i]:
                is_hit = "✓ MATCH" if term in true_terms else ""
                print(f"      - {term}: {prob:.3f}  {is_hit}")
        else:
            print("      (No term exceeded threshold 0.1)")

    # 9. Simple Micro/Macro Precision & Recall at threshold 0.5
    y_test_dense = y_test.toarray()
    tp = np.sum((test_preds_binary == 1) & (y_test_dense == 1))
    fp = np.sum((test_preds_binary == 1) & (y_test_dense == 0))
    fn = np.sum((test_preds_binary == 0) & (y_test_dense == 1))
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    print("-" * 70)
    print(f"QUICK TEST-SET METRICS (at fixed threshold 0.5):")
    print(f"  • True Positives:  {tp}")
    print(f"  • False Positives: {fp}")
    print(f"  • False Negatives: {fn}")
    print(f"  • Micro-Precision: {precision:.3f}")
    print(f"  • Micro-Recall:    {recall:.3f}")
    print(f"  • Micro-F1:        {f1:.3f}")
    print("=" * 70)
    print("Run completed successfully!")


if __name__ == "__main__":
    main()
