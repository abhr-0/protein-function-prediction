#!/usr/bin/env python3
"""Download external data files and record provenance.

Examples:
    # UniProt flat file (GO annotations) + current GO ontology
    .venv/bin/python scripts/download_data.py

    # Also fetch/replace the Swiss-Prot FASTA
    .venv/bin/python scripts/download_data.py --include fasta

    # Re-download everything (moves to a newer release if one exists)
    .venv/bin/python scripts/download_data.py --include dat fasta obo --force
"""

import argparse
import sys
from pathlib import Path

# Ensure project root is in pythonpath when run as a script.
project_root = Path(__file__).resolve().parents[1]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.utils.download import (  # noqa: E402
    download_file,
    fetch_uniprot_release,
    load_manifest,
    save_manifest,
    update_manifest,
    validate_gzip,
)
from src.utils.paths import (  # noqa: E402
    get_data_manifest_path,
    get_go_obo_path,
    get_uniprot_dat_path,
    get_uniprot_path,
)

UNIPROT_BASE = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"
)
GO_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"

SOURCES = {
    "dat": ("uniprot_sprot.dat.gz", f"{UNIPROT_BASE}/uniprot_sprot.dat.gz", get_uniprot_dat_path),
    "fasta": ("uniprot_sprot.fasta.gz", f"{UNIPROT_BASE}/uniprot_sprot.fasta.gz", get_uniprot_path),
    "obo": ("go-basic.obo", GO_OBO_URL, get_go_obo_path),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--include",
        nargs="+",
        choices=sorted(SOURCES),
        default=["dat", "obo"],
        help="Files to download (default: dat obo; fasta is assumed to exist).",
    )
    parser.add_argument(
        "--release",
        default=None,
        help="Override the UniProt release label recorded in the manifest "
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


def main() -> None:
    args = parse_args()

    if args.release:
        release = args.release
    else:
        try:
            release = fetch_uniprot_release() or "unknown"
        except Exception as exc:  # network or certificate failure
            print(f"Warning: could not detect UniProt release: {exc}", file=sys.stderr)
            release = "unknown"

    print(f"UniProt release: {release}")
    manifest_path = get_data_manifest_path()
    manifest = load_manifest(manifest_path)

    for key in args.include:
        name, url, path_getter = SOURCES[key]
        destination = path_getter()

        if args.dry_run:
            action = "would re-download" if args.force else "would fetch/reuse"
            print(f"{action} {name} -> {destination}")
            continue

        print(f"Fetching {name} -> {destination}")
        record = download_file(
            url=url,
            destination=destination,
            name=name,
            release=release,
            force=args.force,
            progress=not args.no_progress,
        )

        if name.endswith(".gz"):
            validate_gzip(destination)

        update_manifest(manifest, record)
        save_manifest(manifest, manifest_path)
        print(
            f"  {record.status}: {record.bytes:,} bytes | sha256={record.sha256[:16]}..."
        )

    if not args.dry_run:
        print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
