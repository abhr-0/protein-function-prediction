"""Tests for project path helpers."""

import pytest

from src.utils.paths import (
    get_annotations_cache_path,
    get_cafa_dir,
    get_cafa_eval_dir,
    get_cafa_eval_obo_path,
    get_cafa_groundtruth_path,
    get_cafa_terms_of_interest_path,
    get_cafa_train_sequences_path,
    get_data_manifest_path,
    get_go_obo_path,
    get_ia_path,
    get_train_terms_path,
    get_uniprot_dat_path,
    get_uniprot_path,
)


def test_uniprot_data_paths():
    """Verify stable filenames for the sequence and annotation sources."""
    assert get_uniprot_path().name == "uniprot_sprot.fasta.gz"
    assert get_uniprot_dat_path().name == "uniprot_sprot.dat.gz"
    assert get_go_obo_path().name == "go-basic.obo"
    assert get_uniprot_dat_path().parent == get_uniprot_path().parent


def test_data_files_are_grouped_by_source():
    """UniProt files share a directory; the GO ontology has its own.

    The ontology is published by the GO Consortium on its own release cadence,
    so co-locating it with UniProt data would imply the wrong provenance.
    """
    assert get_uniprot_path().parent.name == "uniprot"
    assert get_go_obo_path().parent.name == "go"
    assert get_go_obo_path().parent != get_uniprot_path().parent
    assert get_uniprot_path().parent.parent == get_go_obo_path().parent.parent


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


def test_annotations_cache_path_can_include_source_fingerprint():
    """Different source files must not share one cache file."""
    plain = get_annotations_cache_path("2026_04", "cafa")
    keyed = get_annotations_cache_path("2026_04", "cafa", source_fingerprint="abc123")
    other = get_annotations_cache_path("2026_04", "cafa", source_fingerprint="def456")

    assert plain.name == "annotations_sprot_2026_04_cafa.tsv.gz"
    assert keyed.name == "annotations_sprot_2026_04_cafa_abc123.tsv.gz"
    assert len({plain, keyed, other}) == 3


def test_cafa6_layout_separates_train_test_and_eval():
    """Kaggle-sourced bundle and figshare eval bundle must not be mixed."""
    assert get_cafa_dir().name == "cafa6"
    assert get_train_terms_path().parent.name == "train"
    assert get_cafa_train_sequences_path().parent.name == "train"
    assert get_ia_path().parent == get_cafa_eval_dir()
    assert get_cafa_eval_obo_path().parent == get_cafa_eval_dir()
    assert get_cafa_terms_of_interest_path().parent == get_cafa_eval_dir()


def test_cafa6_groundtruth_paths_by_regime():
    """NK/LK/PK regimes map to files; anything else must fail loudly."""
    assert get_cafa_groundtruth_path("NK").name == "May2026_groundtruth_NK.tsv"
    assert get_cafa_groundtruth_path("lk").name == "May2026_groundtruth_LK.tsv"
    assert get_cafa_groundtruth_path("PK_known").name == "May2026_groundtruth_PK_known.tsv"

    with pytest.raises(ValueError, match="Unknown CAFA6 regime"):
        get_cafa_groundtruth_path("XX")
