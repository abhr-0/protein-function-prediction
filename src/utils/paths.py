"""Path utilities for locating project data files and artifacts."""

import re
from pathlib import Path


def get_project_root() -> Path:
    """Return the absolute path to the project root directory."""
    # src/utils/paths.py -> parents[2] is project root
    return Path(__file__).resolve().parents[2]


def get_data_dir() -> Path:
    """Return path to data directory, creating it if needed."""
    data_dir = get_project_root() / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir


def get_uniprot_path() -> Path:
    """Return path to the UniProtKB/Swiss-Prot compressed FASTA."""
    return get_data_dir() / "uniprot_sprot.fasta.gz"


def get_uniprot_dat_path() -> Path:
    """Return path to the UniProtKB/Swiss-Prot flat-file archive (GO annotations)."""
    return get_data_dir() / "uniprot_sprot.dat.gz"


def get_go_obo_path() -> Path:
    """Return path to the Gene Ontology OBO graph (`go-basic.obo`)."""
    return get_data_dir() / "go-basic.obo"


def get_data_manifest_path() -> Path:
    """Return path to the provenance manifest for downloaded external data."""
    return get_artifacts_dir() / "data_manifest.json"


def _slug(value: str) -> str:
    """Return a filesystem-safe, lowercase slug for a manifest key."""
    slug = re.sub(r"[^a-z0-9_-]+", "-", value.strip().lower()).strip("-_")
    return slug or "default"


def get_annotations_cache_path(
    release: str = "2026_04",
    evidence_policy: str = "cafa",
) -> Path:
    """Return path to the cached UniProt GO annotation table.

    Args:
        release: UniProt release identifier, e.g. ``"2026_04"``.
        evidence_policy: Evidence filter name, e.g. ``"cafa"`` or ``"full"``.
    """
    filename = f"annotations_sprot_{_slug(release)}_{_slug(evidence_policy)}.tsv.gz"
    return get_artifacts_dir() / filename


def get_cafa_dir() -> Path:
    """Return path to CAFA dataset directory."""
    return get_data_dir() / "cafa_6_protein_data"


def get_train_terms_path() -> Path:
    """Return path to train_terms.tsv."""
    return get_cafa_dir() / "Train" / "train_terms.tsv"


def get_ia_path() -> Path:
    """Return path to IA.tsv (Information Accretion weights)."""
    return get_cafa_dir() / "IA.tsv"


def get_artifacts_dir() -> Path:
    """Return path to artifacts directory, creating it if needed."""
    artifacts = get_data_dir() / "artifacts"
    artifacts.mkdir(parents=True, exist_ok=True)
    return artifacts
