"""Tests for download and provenance helpers (offline, using file:// URLs)."""

import hashlib
import json
from pathlib import Path

import pytest

from src.utils.download import (
    DownloadRecord,
    download_file,
    file_fingerprint,
    load_manifest,
    manifest_release,
    parse_release_version,
    save_manifest,
    sha256_file,
    update_manifest,
    validate_gzip,
)


def test_parse_release_version():
    """Release notes should yield the YYYY_NN identifier."""
    text = "UniProt Release 2026_04\n\nUniProtKB Release 2026_04 consists of ..."
    assert parse_release_version(text) == "2026_04"
    assert parse_release_version("no release here") is None


def test_sha256_file(tmp_path: Path):
    """Checksums must match a direct hashlib computation."""
    path = tmp_path / "payload.bin"
    payload = b"protein-function-prediction"
    path.write_bytes(payload)

    assert sha256_file(path) == hashlib.sha256(payload).hexdigest()


def test_download_file_copies_and_records(tmp_path: Path):
    """A file:// source should be streamed to the destination."""
    source = tmp_path / "source.bin"
    source.write_bytes(b"ACDEFGHIKLMNPQRSTVWY" * 1000)
    destination = tmp_path / "nested" / "target.bin"

    record = download_file(
        source.as_uri(),
        destination,
        name="target.bin",
        release="2026_04",
        progress=False,
    )

    assert record.status == "downloaded"
    assert record.bytes == source.stat().st_size
    assert record.sha256 == sha256_file(destination)
    assert destination.read_bytes() == source.read_bytes()
    assert record.release == "2026_04"


def test_download_file_reuses_existing_without_force(tmp_path: Path):
    """An existing destination should be checksummed, not overwritten."""
    source = tmp_path / "source.bin"
    source.write_bytes(b"first version")
    destination = tmp_path / "target.bin"

    first = download_file(source.as_uri(), destination, progress=False)
    source.write_bytes(b"second version")
    second = download_file(source.as_uri(), destination, progress=False)

    assert first.status == "downloaded"
    assert second.status == "existing"
    assert destination.read_bytes() == b"first version"
    assert second.sha256 == first.sha256


def test_download_file_force_updates(tmp_path: Path):
    """force=True should replace an existing file and mark it updated."""
    source = tmp_path / "source.bin"
    source.write_bytes(b"first version")
    destination = tmp_path / "target.bin"
    download_file(source.as_uri(), destination, progress=False)

    source.write_bytes(b"second version")
    record = download_file(source.as_uri(), destination, force=True, progress=False)

    assert record.status == "updated"
    assert destination.read_bytes() == b"second version"


def test_download_file_leaves_no_partial_on_error(tmp_path: Path):
    """A failed transfer must not leave a .part file behind."""
    destination = tmp_path / "target.bin"
    missing = (tmp_path / "does-not-exist.bin").as_uri()

    with pytest.raises(Exception):
        download_file(missing, destination, progress=False)

    assert not destination.exists()
    assert not destination.with_name(destination.name + ".part").exists()


def test_manifest_roundtrip_and_merge(tmp_path: Path):
    """Saving then loading should preserve existing entries and add new ones."""
    manifest_path = tmp_path / "manifest.json"
    record = DownloadRecord(
        name="new.tsv.gz",
        url="file:///tmp/new.tsv.gz",
        path="/tmp/new.tsv.gz",
        bytes=123,
        sha256="abc123",
        release="2026_04",
        downloaded_at="2026-11-16T00:00:00+00:00",
        status="downloaded",
    )

    manifest = {"files": {"old.txt": {"sha256": "old"}}}
    update_manifest(manifest, record)
    save_manifest(manifest, manifest_path)

    with open(manifest_path, encoding="utf-8") as handle:
        raw = json.load(handle)
    assert raw["release"] == "2026_04"
    assert "old.txt" in raw["files"]
    assert raw["files"]["new.tsv.gz"]["sha256"] == "abc123"

    loaded = load_manifest(manifest_path)
    assert loaded["files"]["new.tsv.gz"]["bytes"] == 123


def test_load_manifest_missing_returns_empty(tmp_path: Path):
    """Missing manifests should not raise."""
    manifest = load_manifest(tmp_path / "absent.json")
    assert manifest == {"files": {}}


def test_validate_gzip_accepts_and_rejects(tmp_path: Path):
    """Gzip magic-byte validation should catch partial/corrupt downloads."""
    good = tmp_path / "good.gz"
    good.write_bytes(b"\x1f\x8b\x08\x00")
    validate_gzip(good)

    bad = tmp_path / "bad.gz"
    bad.write_bytes(b"not gzip")
    with pytest.raises(ValueError):
        validate_gzip(bad)


def test_manifest_release_prefers_flat_file_entry():
    """The flat-file entry is the authoritative label release."""
    manifest = {
        "release": "2024_01",
        "files": {
            "uniprot_sprot.fasta.gz": {"release": "2025_06"},
            "uniprot_sprot.dat.gz": {"release": "2026_04"},
        },
    }
    assert manifest_release(manifest) == "2026_04"


def test_manifest_release_falls_back_to_fasta_then_manifest_level():
    """Missing entries should fall back rather than return None prematurely."""
    fasta_only = {"files": {"uniprot_sprot.fasta.gz": {"release": "2025_06"}}}
    assert manifest_release(fasta_only) == "2025_06"

    top_level_only = {"release": "2024_01", "files": {}}
    assert manifest_release(top_level_only) == "2024_01"


def test_manifest_release_ignores_unknown():
    """"unknown" must not be mistaken for a real release identifier."""
    manifest = {"release": "unknown", "files": {"uniprot_sprot.dat.gz": {"release": "unknown"}}}
    assert manifest_release(manifest) is None
    assert manifest_release({}) is None


def test_file_fingerprint_tracks_file_identity(tmp_path: Path):
    """Fingerprints must be stable per file and change when it is rewritten."""
    path = tmp_path / "source.dat"
    path.write_bytes(b"first")

    first = file_fingerprint(path)
    assert first == file_fingerprint(path), "same file must fingerprint identically"

    path.write_bytes(b"second version")
    assert file_fingerprint(path) != first

    with pytest.raises(FileNotFoundError):
        file_fingerprint(tmp_path / "missing.dat")
