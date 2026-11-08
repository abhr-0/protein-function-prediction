"""Path utilities for locating project data files and artifacts."""

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
