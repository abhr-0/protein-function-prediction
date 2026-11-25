"""Command-line interface."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

from pfpred import cafa6
from pfpred.labels import CURATED, filter_annotations, summarize
from pfpred.paths import Layout
from pfpred.uniprot import read_annotations, read_fasta

GROUNDTRUTH_KINDS = ("NK", "LK", "PK")

EXPECTED: dict[str, int] = {
    "human.proteins": 20429,
    "human.residues": 11417710,
    "human.annotations": 243091,
    "human.terms": 17317,
    "human.proteins_with_labels": 19375,
    "human.excluded_rows": 26819,
    "human.curated.annotations": 216272,
    "human.curated.proteins": 18630,
    "human.curated.terms": 16442,
    "human.curated.F": 46519,
    "human.curated.P": 96961,
    "human.curated.C": 72792,
    "cafa6.train.proteins": 82404,
    "cafa6.train.annotations": 537027,
    "cafa6.train.F": 128452,
    "cafa6.train.P": 250805,
    "cafa6.train.C": 157770,
    "cafa6.species": 1381,
    "gt.NK": 749,
    "gt.LK": 485,
    "gt.PK": 6951,
    "gt.union": 8023,
    "gt.LK_int_PK": 162,
    "gt.in_train": 7219,
    "cafa6.test_sequences": 224121,
}


def collect_counts(layout: Layout) -> dict[str, int]:
    """Read the raw files and return the verification counts."""
    counts: dict[str, int] = {}

    proteins = 0
    residues = 0
    for protein in read_fasta(layout.human_fasta):
        proteins += 1
        residues += len(protein.sequence)
    counts["human.proteins"] = proteins
    counts["human.residues"] = residues

    all_stats = summarize(read_annotations(layout.human_annotations))
    curated_stats = summarize(
        filter_annotations(read_annotations(layout.human_annotations), CURATED)
    )
    counts["human.annotations"] = all_stats.annotations
    counts["human.terms"] = all_stats.terms
    counts["human.proteins_with_labels"] = all_stats.proteins
    counts["human.excluded_rows"] = all_stats.annotations - curated_stats.annotations
    counts["human.curated.annotations"] = curated_stats.annotations
    counts["human.curated.proteins"] = curated_stats.proteins
    counts["human.curated.terms"] = curated_stats.terms
    for aspect in ("F", "P", "C"):
        counts[f"human.curated.{aspect}"] = curated_stats.per_aspect.get(aspect, 0)

    train_stats = summarize(cafa6.read_terms(layout.cafa6_train_terms))
    counts["cafa6.train.proteins"] = train_stats.proteins
    counts["cafa6.train.annotations"] = train_stats.annotations
    for aspect in ("F", "P", "C"):
        counts[f"cafa6.train.{aspect}"] = train_stats.per_aspect.get(aspect, 0)

    taxonomy = cafa6.read_taxonomy(layout.cafa6_train_taxonomy)
    counts["cafa6.species"] = len(set(taxonomy.values()))

    groundtruth = {
        kind: cafa6.groundtruth_proteins(layout.cafa6_groundtruth(kind))
        for kind in GROUNDTRUTH_KINDS
    }
    union = set().union(*groundtruth.values())
    counts["gt.NK"] = len(groundtruth["NK"])
    counts["gt.LK"] = len(groundtruth["LK"])
    counts["gt.PK"] = len(groundtruth["PK"])
    counts["gt.union"] = len(union)
    counts["gt.LK_int_PK"] = len(groundtruth["LK"] & groundtruth["PK"])
    counts["gt.in_train"] = len(union & taxonomy.keys())

    counts["cafa6.test_sequences"] = sum(
        1 for _ in cafa6.read_sequences(layout.cafa6_test_sequences)
    )
    return counts


def compare(actual: Mapping[str, int], expected: Mapping[str, int] = EXPECTED) -> list[str]:
    """Return one message per missing, unexpected, or drifting key."""
    problems: list[str] = []
    for key in sorted(expected):
        if key not in actual:
            problems.append(f"{key}: missing")
        elif actual[key] != expected[key]:
            problems.append(f"{key}: expected {expected[key]}, got {actual[key]}")
    for key in sorted(set(actual) - set(expected)):
        problems.append(f"{key}: unexpected")
    return problems


def _print_table(counts: Mapping[str, int]) -> None:
    width = max((len(key) for key in counts), default=0)
    for key in sorted(counts):
        print(f"{key:<{width}}  {counts[key]:>9}")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pfpred")
    subparsers = parser.add_subparsers(dest="command", required=True)

    verify = subparsers.add_parser("verify", help="check raw data against documented counts")
    verify.add_argument(
        "--root",
        type=Path,
        default=None,
        help="project root containing data/ (default: PFP_ROOT or the repository)",
    )
    verify.add_argument("--json", action="store_true", help="print counts as JSON")
    verify.add_argument("--no-check", action="store_true", help="print counts without comparing")

    args = parser.parse_args(argv)
    if args.command != "verify":  # pragma: no cover - argparse rejects unknown commands
        parser.error(f"unknown command {args.command!r}")

    layout = Layout(args.root.expanduser().resolve()) if args.root else Layout.from_env()
    counts = collect_counts(layout)

    if args.json:
        print(json.dumps(counts, indent=2, sort_keys=True))
    else:
        _print_table(counts)

    if args.no_check:
        return 0
    problems = compare(counts)
    for problem in problems:
        print(f"error: {problem}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
