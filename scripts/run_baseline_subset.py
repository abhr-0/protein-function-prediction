#!/usr/bin/env python3
"""Run the baseline protein function prediction pipeline.

The default training labels come from UniProt/Swiss-Prot GO annotations. CAFA6
labels remain available as an explicit comparison mode.

Examples:
    # Default: UniProt labels, random 5k proteins, Molecular Function
    .venv/bin/python scripts/run_baseline_subset.py

    # CAFA6 labels for a controlled comparison
    .venv/bin/python scripts/run_baseline_subset.py --label-source cafa --n-proteins 5000

    # Full label space for a tight dev run
    .venv/bin/python scripts/run_baseline_subset.py --n-proteins 2000 --min-term-count 20 --top-terms 0
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

# Ensure project root is in pythonpath when run as a script.
project_root = Path(__file__).resolve().parents[1]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.data.annotations import load_uniprot_go_annotations  # noqa: E402
from src.data.targets import (  # noqa: E402
    apply_vocabulary,
    select_target_terms,
    vocabulary_coverage,
)
from src.data.uniprot import build_annotated_dataset, load_annotations  # noqa: E402
from src.evaluation.metrics import (  # noqa: E402
    macro_prf,
    micro_prf,
    select_threshold_micro_f1,
)
from src.features.pipeline import ProteinFeaturePipeline  # noqa: E402
from src.models.baseline import (  # noqa: E402
    BaselineProteinClassifier,
    build_target_matrix,
    split_by_entry_id,
)
from src.utils.paths import (  # noqa: E402
    get_annotations_cache_path,
    get_artifacts_dir,
    get_train_terms_path,
    get_uniprot_dat_path,
    get_uniprot_path,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train and evaluate the baseline classifier on Swiss-Prot proteins."
    )
    parser.add_argument(
        "--label-source",
        choices=["uniprot", "cafa"],
        default="uniprot",
        help="Training labels: UniProt DR GO lines (default) or CAFA6 train_terms.tsv.",
    )
    parser.add_argument(
        "--annotations-path",
        type=Path,
        default=None,
        help="Override the annotation file (default depends on --label-source).",
    )
    parser.add_argument(
        "--n-proteins",
        type=int,
        default=5000,
        help="Number of annotated proteins to sample (<=0 means all).",
    )
    parser.add_argument(
        "--top-terms",
        type=int,
        default=50,
        help="Maximum target GO terms from the training split (0 means no cap).",
    )
    parser.add_argument(
        "--min-term-count",
        type=int,
        default=1,
        help="Minimum number of training proteins a GO term must appear in.",
    )
    parser.add_argument(
        "--min-terms",
        type=int,
        default=1,
        help="Drop proteins with fewer than this many selected target terms.",
    )
    parser.add_argument(
        "--aspect",
        type=str,
        default="F",
        help="GO aspect: 'F' (MF), 'P' (BP), or 'C' (CC).",
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.2,
        help="Fraction of proteins for the test set (default: 0.2).",
    )
    parser.add_argument(
        "--val-size",
        type=float,
        default=0.1,
        help="Fraction of proteins for the validation set (default: 0.1).",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Seed for the train/val/test split and model (default: 42).",
    )
    parser.add_argument(
        "--sample-seed",
        type=int,
        default=42,
        help="Seed for random protein sampling; -1 disables sampling (first N in FASTA order).",
    )
    parser.add_argument(
        "--evidence-policy",
        choices=["cafa", "full"],
        default="cafa",
        help="UniProt evidence filter: cafa excludes IEA/ND/NR; full keeps IEA.",
    )
    parser.add_argument(
        "--release",
        default="2026_04",
        help="UniProt release label for the annotation cache (default: 2026_04).",
    )
    parser.add_argument(
        "--cache-path",
        type=Path,
        default=None,
        help="Override the UniProt annotation cache path.",
    )
    parser.add_argument(
        "--rebuild-cache",
        action="store_true",
        help="Rebuild the annotation cache from the flat file.",
    )
    return parser.parse_args()


def load_training_annotations(args: argparse.Namespace):
    """Load the annotation table and return it with a human-readable description."""
    if args.label_source == "cafa":
        path = args.annotations_path or get_train_terms_path()
        frame = load_annotations(terms_path=path, aspect=args.aspect)
        return frame, path, "CAFA6 train_terms"

    path = args.annotations_path or get_uniprot_dat_path()
    if not path.exists():
        raise FileNotFoundError(
            f"UniProt flat file not found: {path}\n"
            "Run: .venv/bin/python scripts/download_data.py"
        )

    cache_path = args.cache_path or get_annotations_cache_path(
        release=args.release,
        evidence_policy=args.evidence_policy,
    )
    frame = load_uniprot_go_annotations(
        path,
        evidence_policy=args.evidence_policy,
        aspects=[args.aspect],
        cache_path=cache_path,
        force=args.rebuild_cache,
    )
    return frame, path, f"UniProt release {args.release} ({args.evidence_policy} evidence)"


def main() -> None:
    args = parse_args()

    if args.min_terms < 1:
        raise ValueError("--min-terms must be at least 1")
    if args.min_term_count < 1:
        raise ValueError("--min-term-count must be at least 1")

    max_proteins = args.n_proteins if args.n_proteins > 0 else None
    sample_seed = None if args.sample_seed < 0 else args.sample_seed
    dataset_label = f"{args.n_proteins} proteins" if max_proteins else "FULL DATASET"

    print("=" * 78)
    print("PROTEIN FUNCTION PREDICTION - BASELINE RUN")
    print("=" * 78)
    print(f"Label source:   {args.label_source}")
    print(f"Dataset target: {dataset_label}")
    print(f"Aspect:         {args.aspect}")
    print(f"Test/val split: {args.test_size:.0%} / {args.val_size:.0%}")
    print(f"Sampling seed:  {args.sample_seed}")
    print("-" * 78)

    # 1. Load labels.
    t0 = time.time()
    print("1. Loading annotations...")
    annotations_df, annotations_path, label_description = load_training_annotations(args)
    labeled_ids = set(annotations_df["EntryID"].astype(str).unique())
    print(f"   Label file:   {annotations_path}")
    print(f"   Provenance:   {label_description}")
    print(f"   Annotations:  {len(annotations_df):,} rows across {len(labeled_ids):,} proteins")
    print(f"   Load time:    {time.time() - t0:.2f}s")

    # 2. Build the sequence + all-terms dataset.
    t0 = time.time()
    print("2. Aligning sequences with annotations...")
    df = build_annotated_dataset(
        fasta_path=get_uniprot_path(),
        annotations=annotations_df,
        aspect=args.aspect,
        max_proteins=max_proteins,
        sample_seed=sample_seed,
    )
    matched_ids = set(df["EntryID"].astype(str))
    unmatched_ids = labeled_ids - matched_ids
    if sample_seed is not None and max_proteins is not None:
        print(
            f"   Sampled {len(df):,} of {len(labeled_ids):,} labeled proteins "
            f"({len(unmatched_ids):,} not sampled, or absent from the FASTA)"
        )
    else:
        print(
            f"   Matched {len(df):,} proteins "
            f"({len(unmatched_ids):,} labeled proteins absent from the FASTA)"
        )
    print(f"   Alignment time: {time.time() - t0:.2f}s")

    if len(df) < 10:
        raise ValueError(
            f"Only {len(df)} proteins matched. Check the annotation file and FASTA paths."
        )

    # 3. Split before choosing the label vocabulary.
    print("3. Splitting proteins by EntryID...")
    if args.val_size > 0:
        train_df, val_df, test_df = split_by_entry_id(
            df,
            test_size=args.test_size,
            val_size=args.val_size,
            random_state=args.random_state,
        )
    else:
        train_df, test_df = split_by_entry_id(
            df,
            test_size=args.test_size,
            random_state=args.random_state,
        )
        val_df = df.iloc[0:0].copy()
    print(f"   Train: {len(train_df):,} | Val: {len(val_df):,} | Test: {len(test_df):,}")

    # 4. Choose target terms on the training split only.
    top_k = args.top_terms if args.top_terms > 0 else None
    target_terms = select_target_terms(
        train_df["terms"].tolist(),
        min_count=args.min_term_count,
        top_k=top_k,
    )
    if not target_terms:
        raise ValueError(
            "No target terms survived frequency filtering. "
            "Lower --min-term-count or increase --n-proteins."
        )

    coverage = {
        "train": vocabulary_coverage(train_df["terms"].tolist(), target_terms),
        "val": vocabulary_coverage(val_df["terms"].tolist(), target_terms),
        "test": vocabulary_coverage(test_df["terms"].tolist(), target_terms),
    }
    print(
        f"4. Target vocabulary: {len(target_terms):,} terms "
        f"(train-only, min_count={args.min_term_count}, top_k={top_k})"
    )
    print(
        f"   Test annotation coverage: {coverage['test']['annotation_coverage']:.1%} "
        f"({coverage['test']['proteins_with_target']:,}/{coverage['test']['proteins']:,} proteins)"
    )

    # 5. Restrict terms and drop zero-positive proteins.
    #    Vocabulary was chosen on raw train terms above; apply it to every split
    #    and never re-select it afterwards.
    train_df, train_vocab_stats = apply_vocabulary(
        train_df, target_terms, min_terms=args.min_terms
    )
    val_vocab_stats = None
    if len(val_df):
        val_df, val_vocab_stats = apply_vocabulary(
            val_df, target_terms, min_terms=args.min_terms
        )
    test_df, test_vocab_stats = apply_vocabulary(
        test_df, target_terms, min_terms=args.min_terms
    )
    vocab_stats = {"train": train_vocab_stats, "test": test_vocab_stats}
    if val_vocab_stats is not None:
        vocab_stats["val"] = val_vocab_stats
    print(
        f"5. After term filtering: Train {len(train_df):,} | "
        f"Val {len(val_df):,} | Test {len(test_df):,}"
    )
    print(
        f"   Dropped (fewer than {args.min_terms} in-vocabulary term(s)): "
        f"Train {train_vocab_stats['proteins_dropped']:,} | "
        f"Val {0 if val_vocab_stats is None else val_vocab_stats['proteins_dropped']:,} | "
        f"Test {test_vocab_stats['proteins_dropped']:,}"
    )
    if len(train_df) == 0 or len(test_df) == 0:
        raise ValueError("No proteins remain after term filtering; loosen the vocabulary filters.")

    # 6. Extract features on train only.
    t0 = time.time()
    print("6. Extracting features (AAC + physicochemical + 2-mers)...")
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
    X_val = (
        feature_pipeline.transform(val_df["sequence"].tolist())
        if len(val_df)
        else None
    )
    X_test = feature_pipeline.transform(test_df["sequence"].tolist())
    print(
        f"   Train {X_train.shape} | Val {None if X_val is None else X_val.shape} | "
        f"Test {X_test.shape} in {time.time() - t0:.2f}s"
    )

    # 7. Build targets.
    y_train = build_target_matrix(train_df["terms"].tolist(), target_terms, return_sparse=True)
    y_val = (
        build_target_matrix(val_df["terms"].tolist(), target_terms, return_sparse=True)
        if len(val_df)
        else None
    )
    y_test = build_target_matrix(test_df["terms"].tolist(), target_terms, return_sparse=True)
    print(
        f"7. Targets: train {y_train.shape} "
        f"(sparsity {100 * (1 - y_train.nnz / (y_train.shape[0] * y_train.shape[1])):.1f}%)"
    )

    # 8. Train.
    t0 = time.time()
    print("8. Training BaselineProteinClassifier (one-vs-rest logistic regression)...")
    model = BaselineProteinClassifier(
        target_terms=target_terms,
        n_jobs=-1,
        random_state=args.random_state,
    )
    model.fit(X_train, y_train)
    print(f"   Training time: {time.time() - t0:.2f}s")

    # 9. Tune the threshold on validation, report on test.
    if X_val is not None and y_val is not None and len(val_df):
        val_probs = model.predict_proba(X_val)
        threshold, val_f1 = select_threshold_micro_f1(y_val, val_probs)
        threshold_source = "validation micro-F1"
    else:
        threshold, val_f1 = 0.5, 0.0
        threshold_source = "fixed default"

    test_probs = model.predict_proba(X_test)
    test_pred = (test_probs >= threshold).astype(np.int8)
    micro = micro_prf(y_test, test_pred)
    macro = macro_prf(y_test, test_pred)

    print(f"9. Threshold: {threshold:.2f} ({threshold_source}; val micro-F1={val_f1:.3f})")
    print(
        f"   Test micro: P={micro['precision']:.3f} R={micro['recall']:.3f} "
        f"F1={micro['f1']:.3f} (TP={micro['tp']:,} FP={micro['fp']:,} FN={micro['fn']:,})"
    )
    print(
        f"   Test macro: P={macro['precision']:.3f} R={macro['recall']:.3f} "
        f"F1={macro['f1']:.3f} over {macro['terms']:,} terms"
    )

    # 10. Save model + metrics with provenance.
    tag = f"{args.n_proteins}" if max_proteins else "full"
    artifact_base = f"baseline_model_{args.label_source}_{tag}_{args.aspect}"
    model_path = model.save(get_artifacts_dir() / f"{artifact_base}.joblib")

    metrics = {
        "config": {
            "label_source": args.label_source,
            "label_file": str(annotations_path),
            "label_description": label_description,
            "evidence_policy": args.evidence_policy if args.label_source == "uniprot" else None,
            "release": args.release if args.label_source == "uniprot" else None,
            "aspect": args.aspect,
            "n_proteins_requested": args.n_proteins,
            "sample_seed": args.sample_seed,
            "test_size": args.test_size,
            "val_size": args.val_size,
            "random_state": args.random_state,
            "min_terms": args.min_terms,
            "min_term_count": args.min_term_count,
            "top_terms": args.top_terms,
        },
        "counts": {
            "labeled_proteins": len(labeled_ids),
            "matched_proteins": len(df),
            "unmatched_proteins": len(unmatched_ids),
            "train_proteins": len(train_df),
            "val_proteins": len(val_df),
            "test_proteins": len(test_df),
            "target_terms": len(target_terms),
        },
        "coverage": coverage,
        "vocabulary_filtering": vocab_stats,
        "threshold": {
            "value": threshold,
            "source": threshold_source,
            "validation_micro_f1": val_f1,
        },
        "test": {
            "micro": micro,
            "macro": macro,
        },
    }
    metrics_path = get_artifacts_dir() / f"baseline_metrics_{args.label_source}_{tag}_{args.aspect}.json"
    with open(metrics_path, "w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"10. Saved model:   {model_path}")
    print(f"    Saved metrics: {metrics_path}")

    # 11. Inspection.
    print("-" * 78)
    print("SAMPLE PREDICTIONS (first 3 test proteins):")
    print("-" * 78)
    top_k_preds = model.predict_top_k(X_test[:3], k=5, threshold=0.1)
    for i in range(min(3, len(test_df))):
        entry_id = test_df.iloc[i]["EntryID"]
        name = test_df.iloc[i]["entry_name"]
        org = test_df.iloc[i]["organism"]
        true_terms = test_df.iloc[i]["terms"]
        print(
            f"\nProtein #{i + 1}: {entry_id} ({name}) | "
            f"Length: {test_df.iloc[i]['length']} aa | Organism: {org}"
        )
        print(f"  Ground truth ({len(true_terms)}): {', '.join(true_terms[:6])}")
        print("  Top predictions:")
        for term, prob in top_k_preds[i]:
            marker = " MATCH" if term in true_terms else ""
            print(f"    - {term}: {prob:.3f}{marker}")

    print("=" * 78)
    print("Run completed successfully.")


if __name__ == "__main__":
    main()
