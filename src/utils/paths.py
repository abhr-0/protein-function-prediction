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


def get_uniprot_dir() -> Path:
    """Return the directory holding UniProt-published data files."""
    path = get_data_dir() / "uniprot"
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_go_dir() -> Path:
    """Return the directory holding Gene Ontology data files."""
    path = get_data_dir() / "go"
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_uniprot_path() -> Path:
    """Return path to the UniProtKB/Swiss-Prot compressed FASTA."""
    return get_uniprot_dir() / "uniprot_sprot.fasta.gz"


def get_uniprot_dat_path() -> Path:
    """Return path to the UniProtKB/Swiss-Prot flat-file archive (GO annotations)."""
    return get_uniprot_dir() / "uniprot_sprot.dat.gz"


def get_go_obo_path() -> Path:
    """Return path to the Gene Ontology OBO graph (`go-basic.obo`).

    Lives under ``data/go`` rather than ``data/uniprot``: the ontology is
    published by the GO Consortium on its own release cadence.
    """
    return get_go_dir() / "go-basic.obo"


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
    source_fingerprint: str | None = None,
) -> Path:
    """Return path to the cached UniProt GO annotation table.

    Args:
        release: UniProt release identifier, e.g. ``"2026_04"``.
        evidence_policy: Evidence filter name, e.g. ``"cafa"`` or ``"full"``.
        source_fingerprint: Optional identifier for the exact source file (see
            :func:`src.utils.download.file_fingerprint`). Without it, two
            different flat files sharing a release label would share a cache.
    """
    filename = f"annotations_sprot_{_slug(release)}_{_slug(evidence_policy)}"
    if source_fingerprint:
        filename += f"_{_slug(source_fingerprint)}"
    return get_artifacts_dir() / f"{filename}.tsv.gz"


def get_cafa_dir() -> Path:
    """Return the CAFA6 data directory."""
    return get_data_dir() / "cafa6"


def get_cafa_train_dir() -> Path:
    """Return the CAFA6 training bundle (Kaggle competition files)."""
    return get_cafa_dir() / "train"


def get_cafa_test_dir() -> Path:
    """Return the CAFA6 test superset (Kaggle competition files)."""
    return get_cafa_dir() / "test"


def get_cafa_eval_dir() -> Path:
    """Return the CAFA6 official evaluation bundle (figshare, CC BY 4.0)."""
    return get_cafa_dir() / "eval"


def get_train_terms_path() -> Path:
    """Return path to the CAFA6 ``train_terms.tsv`` labels."""
    return get_cafa_train_dir() / "train_terms.tsv"


def get_cafa_train_sequences_path() -> Path:
    """Return path to CAFA6 ``train_sequences.fasta``.

    Use this for CAFA-mode runs so sequences and labels come from one snapshot.
    """
    return get_cafa_train_dir() / "train_sequences.fasta"


def get_cafa_test_superset_path() -> Path:
    """Return path to CAFA6 ``testsuperset.fasta`` (submission targets)."""
    return get_cafa_test_dir() / "testsuperset.fasta"


def get_cafa_sample_submission_path() -> Path:
    """Return path to CAFA6 ``sample_submission.tsv`` (submission format)."""
    return get_cafa_dir() / "sample_submission.tsv"


def get_ia_path() -> Path:
    """Return path to the official ``IA.tsv`` (Information Accretion weights)."""
    return get_cafa_eval_dir() / "IA.tsv"


def get_cafa_eval_obo_path() -> Path:
    """Return path to CAFA6's own ``go-basic.obo``.

    Use this graph for CAFA6 scoring; use ``get_go_obo_path()`` for our own
    labels. The two are different vintages.
    """
    return get_cafa_eval_dir() / "go-basic.obo"


def get_cafa_groundtruth_path(regime: str = "NK") -> Path:
    """Return path to a CAFA6 ``May2026_groundtruth_<regime>.tsv`` file.

    Args:
        regime: ``"NK"`` (no knowledge), ``"LK"`` (limited knowledge), ``"PK"``
            (partial knowledge), or ``"PK_known"`` (previously public terms).
    """
    regimes = {"NK", "LK", "PK", "PK_known"}
    key = regime.strip().upper() if regime.upper() != "PK_KNOWN" else "PK_known"
    if key not in regimes:
        raise ValueError(f"Unknown CAFA6 regime '{regime}'. Choose from {sorted(regimes)}.")
    return get_cafa_eval_dir() / f"May2026_groundtruth_{key}.tsv"


def get_cafa_terms_of_interest_path() -> Path:
    """Return path to CAFA6's evaluated terms of interest."""
    return get_cafa_eval_dir() / "May2026_groundtruth_terms_of_interest.txt"


def get_artifacts_dir() -> Path:
    """Return path to artifacts directory, creating it if needed."""
    artifacts = get_data_dir() / "artifacts"
    artifacts.mkdir(parents=True, exist_ok=True)
    return artifacts
