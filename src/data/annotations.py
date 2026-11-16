"""UniProt GO annotation parsing, evidence filtering, and caching.

UniProt flat files (``.dat``/``.dat.gz``) carry GO annotations on ``DR GO`` lines:

    DR   GO; GO:0004842; F:ubiquitin-protein transferase activity; IDA:UniProtKB.

This module streams those lines without loading the full archive into memory and
returns a tabular label table compatible with the rest of the pipeline:
``EntryID``, ``term``, ``aspect``, ``evidence``.
"""

import gzip
import os
import tempfile
from collections.abc import Collection, Iterator
from dataclasses import dataclass, replace
from pathlib import Path

import pandas as pd

from src.data.uniprot import normalize_aspect

# Evidence codes that should not be used as training labels.
# CAFA's train_terms appears to exclude IEA; ND/NR carry no biological content.
EVIDENCE_POLICIES: dict[str, frozenset[str]] = {
    "cafa": frozenset({"IEA", "ND", "NR"}),
    "full": frozenset({"ND", "NR"}),
}

ANNOTATION_COLUMNS = ["EntryID", "term", "aspect", "evidence"]


@dataclass(frozen=True)
class GoAnnotation:
    """One direct GO annotation from a UniProt flat-file entry."""

    entry_id: str
    term: str
    aspect: str  # "F", "P", or "C"
    evidence: str  # e.g. "IDA", "IEA", "IMP"


def resolve_excluded_evidence(
    evidence_policy: str | Collection[str] = "cafa",
) -> frozenset[str]:
    """Return the set of evidence codes to exclude for a named policy.

    Args:
        evidence_policy: ``"cafa"`` (exclude IEA/ND/NR), ``"full"`` (exclude
            ND/NR), or an explicit collection of evidence codes to exclude.
    """
    if isinstance(evidence_policy, str):
        key = evidence_policy.strip().lower()
        if key not in EVIDENCE_POLICIES:
            valid = ", ".join(sorted(EVIDENCE_POLICIES))
            raise ValueError(
                f"Unknown evidence policy '{evidence_policy}'. Choose from {valid} "
                "or pass an explicit collection of evidence codes."
            )
        return EVIDENCE_POLICIES[key]
    return frozenset(str(code).strip().upper() for code in evidence_policy)


def parse_uniprot_go_line(line: str) -> GoAnnotation | None:
    """Parse one ``DR   GO;`` line into term, aspect, and evidence code.

    Returns ``None`` for malformed or non-GO lines. The entry accession is not
    part of this line and must be supplied by the caller.
    """
    if not line.startswith("DR   GO;"):
        return None

    body = line[len("DR   ") :].strip()
    try:
        prefix, evidence_field = body.rsplit(";", 1)
        _, term_field, aspect_field = prefix.split(";", 2)
    except ValueError:
        return None

    term = term_field.strip().upper()
    aspect_field = aspect_field.strip()
    evidence = evidence_field.strip().rstrip(".").split(":", 1)[0].strip().upper()

    if not term.startswith("GO:") or not aspect_field or not evidence:
        return None
    aspect = aspect_field[0].upper()
    if aspect not in {"F", "P", "C"}:
        return None

    return GoAnnotation(entry_id="", term=term, aspect=aspect, evidence=evidence)


def stream_uniprot_go_annotations(
    dat_path: str | Path,
    *,
    excluded_evidence: Collection[str] = EVIDENCE_POLICIES["cafa"],
    aspects: Collection[str] | None = None,
) -> Iterator[GoAnnotation]:
    """Stream GO annotations from a UniProt flat file.

    Args:
        dat_path: Path to a ``.dat`` or ``.dat.gz`` file.
        excluded_evidence: Evidence codes to drop.
        aspects: Optional aspects to keep (``"F"``, ``"P"``, ``"C"`` or their
            long forms). ``None`` keeps all aspects.

    Yields:
        GoAnnotation records with the primary accession attached.
    """
    path = Path(dat_path)
    if not path.exists():
        raise FileNotFoundError(f"UniProt flat file not found: {path}")

    excluded = {str(code).strip().upper() for code in excluded_evidence}
    aspect_set = {normalize_aspect(a) for a in aspects} if aspects is not None else None
    opener = gzip.open if path.suffix == ".gz" else open

    current_entry: str | None = None
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("ID   "):
                current_entry = None
            elif line.startswith("AC   "):
                # The AC field wraps across several lines for merged entries,
                # and only the first value of the first line is the primary
                # accession. Later lines carry secondary accessions, which do
                # not appear in the FASTA and would silently drop annotations.
                if current_entry is None:
                    current_entry = line[len("AC   ") :].split(";", 1)[0].strip()
            elif line.startswith("DR   GO;") and current_entry:
                annotation = parse_uniprot_go_line(line)
                if annotation is None:
                    continue
                if annotation.evidence in excluded:
                    continue
                if aspect_set is not None and annotation.aspect not in aspect_set:
                    continue
                yield replace(annotation, entry_id=current_entry)
            elif line.startswith("//"):
                current_entry = None


def load_uniprot_go_annotations(
    dat_path: str | Path,
    *,
    evidence_policy: str | Collection[str] = "cafa",
    aspects: Collection[str] | None = None,
    cache_path: str | Path | None = None,
    force: bool = False,
) -> pd.DataFrame:
    """Load UniProt GO annotations into a tidy DataFrame.

    The cache stores **all aspects** for a given evidence policy; aspect
    filtering happens on read so the same cache can serve MF/BP/CC runs.

    Args:
        dat_path: Path to a ``.dat`` or ``.dat.gz`` file.
        evidence_policy: Named policy (``"cafa"``, ``"full"``) or an explicit
            collection of evidence codes to exclude.
        aspects: Optional aspect filter for the returned frame.
        cache_path: Optional path for a compressed TSV cache.
        force: If true, ignore an existing cache and rebuild it.

    Returns:
        DataFrame with columns ``EntryID``, ``term``, ``aspect``, ``evidence``.
    """
    path = Path(dat_path)
    if not path.exists():
        raise FileNotFoundError(f"UniProt flat file not found: {path}")

    excluded = resolve_excluded_evidence(evidence_policy)
    cache = Path(cache_path) if cache_path is not None else None

    if cache is not None and cache.is_file() and not force:
        frame = _read_annotation_tsv(cache)
    else:
        frame = _build_annotation_frame(
            path,
            excluded_evidence=excluded,
            cache_path=cache,
        )

    if aspects is not None:
        aspect_set = {normalize_aspect(a) for a in aspects}
        frame = frame[frame["aspect"].isin(aspect_set)].reset_index(drop=True)

    return frame


def _build_annotation_frame(
    dat_path: Path,
    *,
    excluded_evidence: Collection[str],
    cache_path: Path | None,
) -> pd.DataFrame:
    """Stream annotations to a temporary TSV, then load it with categoricals."""
    temp_handle = tempfile.NamedTemporaryFile(
        "wt",
        suffix=".tsv",
        delete=False,
        encoding="utf-8",
    )
    temp_path = Path(temp_handle.name)
    try:
        with temp_handle as handle:
            handle.write("\t".join(ANNOTATION_COLUMNS) + "\n")
            for annotation in stream_uniprot_go_annotations(
                dat_path,
                excluded_evidence=excluded_evidence,
            ):
                handle.write(
                    f"{annotation.entry_id}\t{annotation.term}\t"
                    f"{annotation.aspect}\t{annotation.evidence}\n"
                )

        frame = _read_annotation_tsv(temp_path)
        frame = frame.drop_duplicates().reset_index(drop=True)

        if cache_path is not None:
            _save_annotation_cache(frame, cache_path)
        return frame
    finally:
        if temp_path.exists():
            temp_path.unlink()


def _read_annotation_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        dtype={
            "EntryID": "category",
            "term": "category",
            "aspect": "category",
            "evidence": "category",
        },
    )


def _save_annotation_cache(frame: pd.DataFrame, cache_path: Path) -> Path:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cache_path.with_name(cache_path.name + ".tmp")
    compression = "gzip" if cache_path.suffix == ".gz" else None
    frame.to_csv(tmp_path, sep="\t", index=False, compression=compression)
    os.replace(tmp_path, cache_path)
    return cache_path


def summarize_annotations(frame: pd.DataFrame) -> dict:
    """Summarize annotation coverage, aspect balance, and evidence codes.

    Args:
        frame: Tidy annotation frame with ``EntryID``, ``term``, ``aspect``,
            ``evidence`` columns.

    Returns:
        JSON-serializable dictionary with totals and per-aspect statistics.
    """
    summary: dict = {
        "annotations": int(len(frame)),
        "proteins": int(frame["EntryID"].nunique()) if len(frame) else 0,
        "terms": int(frame["term"].nunique()) if len(frame) else 0,
        "aspects": {},
        "evidence": {},
    }

    if len(frame) == 0:
        return summary

    for aspect, group in frame.groupby("aspect", observed=True):
        per_protein = group.groupby("EntryID", observed=True).size()
        summary["aspects"][str(aspect)] = {
            "annotations": int(len(group)),
            "proteins": int(group["EntryID"].nunique()),
            "terms": int(group["term"].nunique()),
            "terms_per_protein_mean": round(float(per_protein.mean()), 3),
            "terms_per_protein_median": float(per_protein.median()),
            "terms_per_protein_max": int(per_protein.max()),
        }

    summary["evidence"] = {
        str(code): int(count)
        for code, count in frame["evidence"].value_counts().items()
    }
    return summary
