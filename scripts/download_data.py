#!/usr/bin/env python3
"""Download external data files and record provenance.

Examples:
    # Default: Swiss-Prot labels + sequences + GO ontology
    .venv/bin/python scripts/download_data.py

    # Add the official CAFA6 evaluation bundle (figshare, CC BY 4.0)
    .venv/bin/python scripts/download_data.py --include cafa6

    # Everything
    .venv/bin/python scripts/download_data.py --include all

    # Re-download (moves to a newer release if one exists)
    .venv/bin/python scripts/download_data.py --include dat fasta obo --force
"""

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

# Ensure project root is in pythonpath when run as a script.
project_root = Path(__file__).resolve().parents[1]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.utils.download import (  # noqa: E402
    download_file,
    fetch_uniprot_release,
    go_obo_release,
    load_manifest,
    save_manifest,
    update_manifest,
    validate_gzip,
)
from src.utils.paths import (  # noqa: E402
    get_cafa_dir,
    get_cafa_eval_dir,
    get_data_manifest_path,
    get_go_obo_path,
    get_uniprot_dat_path,
    get_uniprot_path,
)

UNIPROT_BASE = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"
)
GO_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"

# Official CAFA6 evaluation data, published by the competition organizers on
# figshare under CC BY 4.0 (DOI 10.6084/m9.figshare.32861969). During the
# competition these files were secret; they were released after it closed, so
# CAFA6 can now be scored locally.
CAFA6_FIGSHARE = "https://ndownloader.figshare.com/files"
CAFA6_RELEASE = "CAFA6 May2026"
CAFA6_EVAL_FILES = [
    # (filename, figshare file id, md5)
    ("IA.tsv", "66194243", "368d7288f0d99dfe6aee594a7bc9d9a6"),
    ("go-basic.obo", "66194237", "54799f704d1dd94882f4ac18ff3cdb7c"),
    ("test_sequences.fasta", "66194249", "20964bca88cd2d39bc5803a352b25c67"),
    ("May2026_groundtruth_NK.tsv", "66194255", "8bb4e48c6afcb457bfece13a5b5f9d88"),
    ("May2026_groundtruth_LK.tsv", "66194258", "a8acf3765083760fdd5d6681cd60dd15"),
    ("May2026_groundtruth_PK.tsv", "66194252", "10afc938b08dd62f6b759ec3ccc01c4d"),
    ("May2026_groundtruth_PK_known.tsv", "66194240", "96deb2067c580e13420c045a6ea1439f"),
    ("May2026_groundtruth_terms_of_interest.txt", "66194246", "84ff532f40e85e719f30af68caa89763"),
    ("README.md", "66194234", "d0bea65f55e9f39a1defe78f904765f4"),
]


@dataclass(frozen=True)
class Resource:
    """One downloadable file and the provenance needed to verify it."""

    name: str
    url: str
    path_getter: Callable[[], Path]
    source: str
    md5: str | None = None
    release: str | None = None  # None means: derive per source


def _cafa6_eval_resources() -> list[Resource]:
    return [
        Resource(
            name=filename,
            url=f"{CAFA6_FIGSHARE}/{file_id}",
            path_getter=lambda filename=filename: get_cafa_eval_dir() / filename,
            source="cafa6-figshare",
            md5=md5,
            release=CAFA6_RELEASE,
        )
        for filename, file_id, md5 in CAFA6_EVAL_FILES
    ]


GROUPS: dict[str, list[Resource]] = {
    # Sequences first: they are needed by everything, and are cheap only
    # relative to the flat file.
    "fasta": [
        Resource(
            "uniprot_sprot.fasta.gz",
            f"{UNIPROT_BASE}/uniprot_sprot.fasta.gz",
            get_uniprot_path,
            source="uniprot",
        )
    ],
    "dat": [
        Resource(
            "uniprot_sprot.dat.gz",
            f"{UNIPROT_BASE}/uniprot_sprot.dat.gz",
            get_uniprot_dat_path,
            source="uniprot",
        )
    ],
    "obo": [Resource("go-basic.obo", GO_OBO_URL, get_go_obo_path, source="go")],
    "cafa6": _cafa6_eval_resources(),
}

DEFAULT_GROUPS = ["dat", "fasta", "obo"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--include",
        nargs="+",
        choices=sorted(GROUPS) + ["all"],
        default=DEFAULT_GROUPS,
        help=f"Groups to download (default: {' '.join(DEFAULT_GROUPS)}).",
    )
    parser.add_argument(
        "--release",
        default=None,
        help="Override the UniProt release label recorded for UniProt files "
        "(default: auto-detect from the FTP release notes).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even if they already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be downloaded without transferring anything.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Suppress progress output.",
    )
    return parser.parse_args()


def _resolve_release(args: argparse.Namespace, resource: Resource) -> str:
    """Return the release label for a resource, per publisher.

    UniProt files take the UniProt release; GO versions itself and is filled in
    from the ontology header after download; CAFA6 declares its own label.
    """
    if resource.release is not None:
        return resource.release
    if resource.source == "go":
        return "unknown"
    return args.release


def main() -> None:
    args = parse_args()

    groups = sorted(GROUPS) if "all" in args.include else args.include

    if args.release:
        uniprot_release = args.release
    elif "dat" in groups or "fasta" in groups:
        try:
            uniprot_release = fetch_uniprot_release() or "unknown"
        except Exception as exc:  # network or certificate failure
            print(f"Warning: could not detect UniProt release: {exc}", file=sys.stderr)
            uniprot_release = "unknown"
    else:
        uniprot_release = "unknown"
    args.release = uniprot_release

    print(f"UniProt release: {uniprot_release}")
    manifest_path = get_data_manifest_path()
    manifest = load_manifest(manifest_path)

    for group in groups:
        for resource in GROUPS[group]:
            destination = resource.path_getter()
            release = _resolve_release(args, resource)

            if args.dry_run:
                print(f"would fetch/reuse {resource.name} -> {destination}")
                continue

            # Publisher-prefixed key: "go-basic.obo" exists in both the GO and
            # CAFA6 bundles, and a shared manifest key would overwrite one.
            key = f"{resource.source}/{resource.name}"
            print(f"Fetching {key} -> {destination}")
            record = download_file(
                url=resource.url,
                destination=destination,
                name=key,
                release=release,
                source=resource.source,
                expected_md5=resource.md5,
                force=args.force,
                progress=not args.no_progress,
            )

            if resource.name.endswith(".gz"):
                validate_gzip(destination)
            if resource.source == "go" and record.release == "unknown":
                # The ontology versions itself via its header date.
                detected = go_obo_release(destination)
                if detected:
                    record = replace_release(record, detected)

            update_manifest(manifest, record)
            save_manifest(manifest, manifest_path)
            print(
                f"  {record.status}: {record.bytes:,} bytes | "
                f"sha256={record.sha256[:16]}... | release={record.release}"
                + (" | md5 verified" if record.verified else "")
            )

    if not args.dry_run:
        print(f"Manifest: {manifest_path}")


def replace_release(record, release: str):
    """Return a copy of ``record`` with a different release label."""
    from dataclasses import replace

    return replace(record, release=release)


if __name__ == "__main__":
    main()
