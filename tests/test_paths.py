"""Tests for project path helpers."""

from src.utils.paths import (
    get_annotations_cache_path,
    get_data_manifest_path,
    get_go_obo_path,
    get_uniprot_dat_path,
    get_uniprot_path,
)


def test_uniprot_data_paths():
    """Verify stable filenames for the sequence and annotation sources."""
    assert get_uniprot_path().name == "uniprot_sprot.fasta.gz"
    assert get_uniprot_dat_path().name == "uniprot_sprot.dat.gz"
    assert get_go_obo_path().name == "go-basic.obo"
    assert get_uniprot_dat_path().parent == get_uniprot_path().parent


def test_data_manifest_path_is_under_artifacts():
    """The provenance manifest should live with other generated artifacts."""
    manifest = get_data_manifest_path()
    assert manifest.name == "data_manifest.json"
    assert manifest.parent.name == "artifacts"


def test_annotations_cache_path_uses_release_and_policy():
    """Cache names should be deterministic and filesystem-safe."""
    path = get_annotations_cache_path("2026_04", "cafa")
    assert path.name == "annotations_sprot_2026_04_cafa.tsv.gz"
    assert path.parent.name == "artifacts"


def test_annotations_cache_path_sanitizes_policy():
    """Spaces and slashes must not leak into the filename."""
    path = get_annotations_cache_path("2026_04", "No IEA/ND")
    assert path.name == "annotations_sprot_2026_04_no-iea-nd.tsv.gz"
