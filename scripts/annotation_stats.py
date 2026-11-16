#!/usr/bin/env python3
"""Report GO annotation coverage and evidence distribution for Swiss-Prot.

Examples:
    .venv/bin/python scripts/annotation_stats.py --aspect F
    .venv/bin/python scripts/annotation_stats.py --evidence-policy full
    .venv/bin/python scripts/annotation_stats.py --json data/artifacts/annotation_stats.json
"""

import argparse
import json
import sys
from pathlib import Path

# Ensure project root is in pythonpath when run as a script.
project_root = Path(__file__).resolve().parents[1]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.data.annotations import (  # noqa: E402
    load_uniprot_go_annotations,
    summarize_annotations,
)
from src.utils.download import file_fingerprint  # noqa: E402
from src.utils.paths import (  # noqa: E402
    get_annotations_cache_path,
    get_uniprot_dat_path,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dat-path",
        type=Path,
        default=get_uniprot_dat_path(),
        help="UniProt flat file (.dat or .dat.gz).",
    )
    parser.add_argument(
        "--release",
        default="2026_04",
        help="Release label used for the annotation cache filename.",
    )
    parser.add_argument(
        "--evidence-policy",
        default="cafa",
        choices=["cafa", "full"],
        help="Evidence filter: cafa excludes IEA/ND/NR; full keeps IEA.",
    )
    parser.add_argument(
        "--aspect",
        default=None,
        choices=["F", "P", "C", "MF", "BP", "CC"],
        help="Restrict the report to one aspect (default: all).",
    )
    parser.add_argument(
        "--cache-path",
        type=Path,
        default=None,
        help="Override the annotation cache path.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rebuild the annotation cache from the flat file.",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=None,
        help="Optional path to write the summary as JSON.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cache_path = args.cache_path or get_annotations_cache_path(
        release=args.release,
        evidence_policy=args.evidence_policy,
        source_fingerprint=file_fingerprint(args.dat_path),
    )

    print(f"Flat file:       {args.dat_path}")
    print(f"Evidence policy: {args.evidence_policy}")
    print(f"Cache:           {cache_path}")
    if args.aspect:
        print(f"Aspect filter:   {args.aspect}")

    frame = load_uniprot_go_annotations(
        args.dat_path,
        evidence_policy=args.evidence_policy,
        aspects=[args.aspect] if args.aspect else None,
        cache_path=cache_path,
        force=args.force,
    )
    summary = summarize_annotations(frame)

    print()
    print(f"Annotations:      {summary['annotations']:,}")
    print(f"Unique proteins:  {summary['proteins']:,}")
    print(f"Unique GO terms:  {summary['terms']:,}")

    if summary["aspects"]:
        print()
        print("By aspect:")
        print(
            f"  {'aspect':<8}{'annotations':>14}{'proteins':>12}{'terms':>10}"
            f"{'terms/protein':>16}"
        )
        for aspect, stats in sorted(summary["aspects"].items()):
            print(
                f"  {aspect:<8}{stats['annotations']:>14,}{stats['proteins']:>12,}"
                f"{stats['terms']:>10,}{stats['terms_per_protein_mean']:>16.2f}"
            )

    if summary["evidence"]:
        print()
        print("Evidence codes:")
        for code, count in sorted(
            summary["evidence"].items(), key=lambda item: item[1], reverse=True
        ):
            print(f"  {code:<8}{count:>12,}")

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
        print(f"\nWrote {args.json}")


if __name__ == "__main__":
    main()
