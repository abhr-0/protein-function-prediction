"""Download helpers and provenance manifest for external data files.

The project consumes large public files (UniProt flat files, GO ontology
snapshots). This module streams downloads to disk, computes checksums, and
records provenance in a JSON manifest so experiments can be reproduced.
"""

import hashlib
import json
import os
import re
import ssl
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

USER_AGENT = "protein-function-prediction/0.1"
CHUNK_SIZE = 1 << 20
PROGRESS_STEP = 100 << 20

UNIPROT_RELEASE_NOTES_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt"
)
_RELEASE_RE = re.compile(r"UniProt Release\s+(\d{4}_\d{2})")


@dataclass(frozen=True)
class DownloadRecord:
    """Provenance record for one downloaded (or reused) file."""

    name: str
    url: str
    path: str
    bytes: int
    sha256: str
    release: str
    downloaded_at: str
    status: str  # "downloaded", "updated", or "existing"
    md5: str | None = None
    source: str | None = None
    verified: str | None = None  # hash algorithm checked against upstream


def utc_now_iso() -> str:
    """Return the current UTC timestamp as an ISO-8601 string."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: str | Path, chunk_size: int = CHUNK_SIZE) -> str:
    """Return the SHA-256 hex digest of a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_fingerprint(path: str | Path) -> str:
    """Return a short fingerprint identifying one local file version.

    Combines the resolved path, size, and modification time. This is cheap
    enough to compute on every pipeline start (no full content hash) and
    changes whenever the file is replaced or regenerated, which is what cache
    keys need.
    """
    resolved = Path(path).resolve()
    stat = resolved.stat()
    digest = hashlib.sha256(
        f"{resolved}|{stat.st_size}|{stat.st_mtime_ns}".encode()
    ).hexdigest()
    return digest[:10]


def parse_release_version(text: str) -> str | None:
    """Extract a UniProt release identifier from release notes text."""
    match = _RELEASE_RE.search(text)
    return match.group(1) if match else None


def fetch_uniprot_release(timeout: int = 30) -> str | None:
    """Fetch the current UniProt release identifier from the FTP server."""
    request = Request(UNIPROT_RELEASE_NOTES_URL, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=timeout, context=_ssl_context()) as response:
        text = response.read().decode("utf-8", errors="replace")
    return parse_release_version(text)


def validate_gzip(path: str | Path) -> None:
    """Raise ValueError if a file does not start with the gzip magic bytes."""
    with open(path, "rb") as handle:
        magic = handle.read(2)
    if magic != b"\x1f\x8b":
        raise ValueError(f"Not a gzip file (bad magic bytes): {path}")


def download_file(
    url: str,
    destination: str | Path,
    *,
    name: str | None = None,
    release: str = "unknown",
    source: str | None = None,
    expected_md5: str | None = None,
    force: bool = False,
    timeout: int = 120,
    progress: bool = True,
) -> DownloadRecord:
    """Download ``url`` to ``destination`` and return a provenance record.

    If the destination already exists and ``force`` is false, the file is
    reused and its checksum is recorded instead of re-downloading it.
    Downloads are written to a temporary ``.part`` file and moved into place
    only after the transfer completes.

    Args:
        expected_md5: Upstream-published md5. When given, the transfer is
            verified and a mismatch deletes the partial file and raises.
        source: Provenance label for the publisher, e.g. ``"uniprot"``,
            ``"go"``, or ``"cafa6-figshare"``. Sources have independent
            release cadences, so this is recorded per file.
    """
    dest = Path(destination)
    dest.parent.mkdir(parents=True, exist_ok=True)
    filename = name or dest.name
    existed_before = dest.is_file() and dest.stat().st_size > 0

    if existed_before and not force:
        data = dest.read_bytes()
        actual_md5 = hashlib.md5(data).hexdigest()
        if expected_md5 is not None and actual_md5 != expected_md5:
            raise ValueError(
                f"md5 mismatch for existing {dest}: expected {expected_md5}, "
                f"got {actual_md5}. Re-run with force=True."
            )
        return DownloadRecord(
            name=filename,
            url=url,
            path=str(dest),
            bytes=dest.stat().st_size,
            sha256=sha256_file(dest),
            release=release,
            downloaded_at=utc_now_iso(),
            status="existing",
            md5=actual_md5,
            source=source,
            verified="md5" if expected_md5 else None,
        )

    part = dest.with_name(dest.name + ".part")
    if part.exists():
        part.unlink()

    digest = hashlib.sha256()
    md5_digest = hashlib.md5()
    size = 0
    next_progress = PROGRESS_STEP
    try:
        with _open_url(url, timeout) as response, open(part, "wb") as handle:
            while True:
                chunk = response.read(CHUNK_SIZE)
                if not chunk:
                    break
                handle.write(chunk)
                digest.update(chunk)
                md5_digest.update(chunk)
                size += len(chunk)
                if progress and size >= next_progress:
                    print(
                        f"  ... {filename}: {size / (1 << 20):,.0f} MiB",
                        file=sys.stderr,
                    )
                    next_progress += PROGRESS_STEP
        actual_md5 = md5_digest.hexdigest()
        if expected_md5 is not None and actual_md5 != expected_md5:
            raise ValueError(
                f"md5 mismatch for {dest}: expected {expected_md5}, got {actual_md5}"
            )
        os.replace(part, dest)
    except Exception:
        if part.exists():
            part.unlink()
        raise

    return DownloadRecord(
        name=filename,
        url=url,
        path=str(dest),
        bytes=size,
        sha256=digest.hexdigest(),
        release=release,
        downloaded_at=utc_now_iso(),
        status="updated" if existed_before else "downloaded",
        md5=actual_md5,
        source=source,
        verified="md5" if expected_md5 else None,
    )


def load_manifest(path: str | Path) -> dict:
    """Load a provenance manifest, returning an empty one if absent."""
    manifest_path = Path(path)
    if not manifest_path.is_file():
        return {"files": {}}

    with open(manifest_path, encoding="utf-8") as handle:
        manifest = json.load(handle)
    manifest.setdefault("files", {})
    return manifest


def update_manifest(manifest: dict, record: DownloadRecord) -> dict:
    """Insert or replace a download record in a manifest in place."""
    manifest.setdefault("files", {})[record.name] = asdict(record)
    manifest["updated_at"] = utc_now_iso()
    if record.release and record.release != "unknown":
        manifest["release"] = record.release
    return manifest


def save_manifest(manifest: dict, path: str | Path) -> Path:
    """Write a manifest atomically and return its path."""
    manifest_path = Path(path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = manifest_path.with_name(manifest_path.name + ".tmp")
    with open(tmp_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp_path, manifest_path)
    return manifest_path


def manifest_release(manifest: dict) -> str | None:
    """Return the UniProt release recorded in a provenance manifest.

    Prefers the entry for the Swiss-Prot flat file, falls back to the FASTA
    entry, then to the manifest-level release. Returns ``None`` when nothing
    usable is recorded, so callers never mistake ``"unknown"`` for a version.
    """
    entries = manifest.get("files", {})
    for wanted in ("uniprot_sprot.dat.gz", "uniprot_sprot.fasta.gz"):
        for key, entry in entries.items():
            if not (key == wanted or key.endswith(f"/{wanted}")):
                continue
            release = entry.get("release")
            if release and release != "unknown":
                return str(release)

    release = manifest.get("release")
    if release and release != "unknown":
        return str(release)
    return None



def _open_url(url: str, timeout: int):
    request = Request(url, headers={"User-Agent": USER_AGENT})
    if urlparse(url).scheme == "https":
        return urlopen(request, timeout=timeout, context=_ssl_context())
    return urlopen(request, timeout=timeout)


def _ssl_context() -> ssl.SSLContext:
    cafile = _find_ca_file()
    return ssl.create_default_context(cafile=cafile) if cafile else ssl.create_default_context()


def _find_ca_file() -> str | None:
    """Locate a CA bundle, working around minimal Nix environments."""
    candidates: list[str | None] = [
        os.environ.get("SSL_CERT_FILE"),
        "/etc/ssl/certs/ca-certificates.crt",
        "/etc/ssl/cert.pem",
        "/etc/ssl/certs/ca-bundle.crt",
    ]
    try:
        import certifi  # type: ignore

        candidates.insert(1, certifi.where())
    except ImportError:
        pass

    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return candidate
    return None


def go_obo_release(path: str | Path, header_lines: int = 40) -> str | None:
    """Return the ontology release date from an OBO header.

    Current GO snapshots carry ``data-version: releases/2026-07-26``; older ones
    carry ``date: 08.08.2026 12:00``. Either way the ontology versions itself,
    independently of any UniProt release. Returns an ISO date or ``None``.
    """
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            for _ in range(header_lines):
                line = handle.readline()
                if not line or line.startswith("["):
                    break
                if line.startswith("data-version:"):
                    raw = line.split(":", 1)[1].strip()
                    candidate = raw.rsplit("/", 1)[-1]
                    year, month, day = candidate.split("-")
                    if len(year) == 4:
                        return f"{year}-{month}-{day}"
                if line.startswith("date:"):
                    raw = line.split(":", 1)[1].strip().split()[0]
                    day, month, year = raw.split(".")
                    return f"{year}-{month}-{day}"
    except (OSError, ValueError):
        return None
    return None
