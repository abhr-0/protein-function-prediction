"""UniProtKB/Swiss-Prot data streaming, header parsing, and annotation joining."""

import gzip
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
from Bio import SeqIO


@dataclass
class UniProtRecord:
    """Represents a parsed UniProt sequence record with structured metadata."""

    entry_id: str  # Accession ID (e.g. 'Q6GZX4')
    entry_name: str  # Identifier (e.g. '001R_FRG3G')
    description: str
    organism: str | None
    taxon_id: int | None
    gene_name: str | None
    protein_evidence: int | None
    sequence_version: int | None
    sequence: str
    length: int


# Standard UniProt header regex pattern
_UNIPROT_HEADER_RE = re.compile(
    r"^(?:sp|tr)\|(?P<accession>[A-Za-z0-9_-]+)\|(?P<entry_name>\S+)\s+"
    r"(?P<description>.*?)"
    r"(?:\s+OS=(?P<organism>.*?))?"
    r"(?:\s+OX=(?P<taxon_id>\d+))?"
    r"(?:\s+GN=(?P<gene_name>\S+))?"
    r"(?:\s+PE=(?P<protein_evidence>\d+))?"
    r"(?:\s+SV=(?P<sequence_version>\d+))?$"
)


def parse_uniprot_header(header: str) -> dict:
    """Parse a UniProt FASTA header into structured fields.

    Handles standard format:
    >sp|Accession|EntryName Description OS=Organism OX=TaxonID GN=Gene PE=Evidence SV=Version
    """
    cleaned = header.lstrip(">").strip()
    match = _UNIPROT_HEADER_RE.match(cleaned)

    if match:
        data = match.groupdict()
        return {
            "entry_id": data["accession"],
            "entry_name": data["entry_name"],
            "description": data["description"].strip(),
            "organism": data["organism"],
            "taxon_id": int(data["taxon_id"]) if data["taxon_id"] else None,
            "gene_name": data["gene_name"],
            "protein_evidence": int(data["protein_evidence"]) if data["protein_evidence"] else None,
            "sequence_version": int(data["sequence_version"]) if data["sequence_version"] else None,
        }

    # Fallback if optional fields don't strictly match regex
    first_token = cleaned.split()[0] if cleaned else ""
    parts = first_token.split("|")
    accession = parts[1] if len(parts) > 1 else first_token
    entry_name = parts[2] if len(parts) > 2 else ""
    description = cleaned[len(first_token) :].strip()

    return {
        "entry_id": accession,
        "entry_name": entry_name,
        "description": description,
        "organism": None,
        "taxon_id": None,
        "gene_name": None,
        "protein_evidence": None,
        "sequence_version": None,
    }


def stream_uniprot_fasta(
    fasta_path: str | Path,
    max_records: int | None = None,
) -> Iterator[UniProtRecord]:
    """Stream UniProt records from a plain or gzip-compressed FASTA file.

    Memory-efficient generator that yields one UniProtRecord at a time without
    loading the full database into RAM.
    """
    path = Path(fasta_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    opener = gzip.open if str(path).endswith(".gz") else open

    with opener(path, "rt", encoding="utf-8") as handle:
        count = 0
        for record in SeqIO.parse(handle, "fasta"):
            parsed = parse_uniprot_header(record.description)
            seq_str = str(record.seq).upper()

            yield UniProtRecord(
                entry_id=parsed["entry_id"],
                entry_name=parsed["entry_name"],
                description=parsed["description"],
                organism=parsed["organism"],
                taxon_id=parsed["taxon_id"],
                gene_name=parsed["gene_name"],
                protein_evidence=parsed["protein_evidence"],
                sequence_version=parsed["sequence_version"],
                sequence=seq_str,
                length=len(seq_str),
            )

            count += 1
            if max_records is not None and count >= max_records:
                break


def normalize_aspect(aspect: str) -> str:
    """Normalize aspect input to standard single-letter code ('F', 'P', 'C')."""
    mapping = {
        "F": "F",
        "MF": "F",
        "MOLECULAR_FUNCTION": "F",
        "P": "P",
        "BP": "P",
        "BIOLOGICAL_PROCESS": "P",
        "C": "C",
        "CC": "C",
        "CELLULAR_COMPONENT": "C",
    }
    key = aspect.strip().upper()
    if key not in mapping:
        raise ValueError(f"Unknown aspect '{aspect}'. Must be one of 'F' (MF), 'P' (BP), or 'C' (CC).")
    return mapping[key]


def load_annotations(
    terms_path: str | Path,
    aspect: str | None = None,
) -> pd.DataFrame:
    """Load GO annotations from a CAFA-style ``train_terms.tsv`` table.

    This is pure I/O: there is deliberately **no** term filtering here. Term
    vocabulary must be chosen on the training split only, via
    :func:`src.data.targets.select_target_terms`, otherwise test proteins leak
    into the label space.

    Args:
        terms_path: Path to a ``train_terms.tsv``-style file with columns
            ``EntryID``, ``term``, ``aspect``.
        aspect: Optional aspect filter ('F'/'MF', 'P'/'BP', 'C'/'CC').

    Returns:
        DataFrame with the loaded annotations.
    """
    terms_file = Path(terms_path)
    if not terms_file.exists():
        raise FileNotFoundError(f"Annotations file not found: {terms_file}")

    df = pd.read_csv(terms_file, sep="\t")

    if aspect is not None:
        normalized = normalize_aspect(aspect)
        df = df[df["aspect"] == normalized]

    return df


def load_information_accretion(ia_path: str | Path) -> dict[str, float]:
    """Load Information Accretion (IA) weights from IA.tsv.

    Returns:
        dict mapping GO term ID (e.g. 'GO:0005515') to IA weight (float).
    """
    file_path = Path(ia_path)
    if not file_path.exists():
        raise FileNotFoundError(f"IA file not found: {file_path}")

    ia_df = pd.read_csv(file_path, sep="\t", header=None, names=["term", "weight"])
    return dict(zip(ia_df["term"], ia_df["weight"]))


def build_annotated_dataset(
    fasta_path: str | Path,
    annotations: pd.DataFrame,
    target_terms: Sequence[str] | None = None,
    *,
    aspect: str | None = None,
    max_proteins: int | None = None,
    sample_seed: int | None = None,
    entry_col: str = "EntryID",
    term_col: str = "term",
    aspect_col: str = "aspect",
) -> pd.DataFrame:
    """Align protein sequences with their GO annotations.

    This builder does alignment only. It applies no vocabulary filtering and no
    zero-positive filtering, because both depend on a target vocabulary that
    must be chosen on the training split *after* the split is made. Use
    :func:`src.data.targets.select_target_terms` to choose the vocabulary and
    :func:`src.data.targets.apply_vocabulary` to restrict the splits.

    Args:
        fasta_path: Path to plain or gzip-compressed FASTA.
        annotations: Tidy annotation frame.
        target_terms: Optional vocabulary to restrict each protein's terms to.
            ``None`` keeps every annotated term.
        aspect: Optional aspect filter applied to ``annotations``.
        max_proteins: Maximum number of proteins to return. ``None`` keeps all.
        sample_seed: If set, randomly sample up to ``max_proteins`` proteins
            across the full labeled set instead of taking the first N records
            in FASTA order.
        entry_col: Accession column name in ``annotations``.
        term_col: GO term column name in ``annotations``.
        aspect_col: Aspect column name in ``annotations``.

    Returns:
        DataFrame with columns ``EntryID``, ``entry_name``, ``organism``,
        ``sequence``, ``length``, ``terms``.
    """
    if aspect is not None:
        normalized = normalize_aspect(aspect)
        annotations = annotations[annotations[aspect_col] == normalized]

    term_set = set(target_terms) if target_terms is not None else None

    # Cast to plain str: annotation sources may store terms as categoricals
    # (the UniProt cache does), and pandas cannot regroup a list-valued
    # aggregation into a categorical column.
    pairs = annotations[[entry_col, term_col]].astype(str)

    terms_grouped: dict[str, list[str]] = {}
    grouped = pairs.groupby(entry_col, observed=True)[term_col].agg(
        lambda values: sorted(set(values))
    )
    for entry_id, terms in grouped.items():
        selected = [term for term in terms if term_set is None or term in term_set]
        if selected:
            terms_grouped[str(entry_id)] = selected

    target_protein_ids = set(terms_grouped)
    if sample_seed is not None and max_proteins is not None:
        candidate_ids = sorted(target_protein_ids)
        if len(candidate_ids) > max_proteins:
            rng = np.random.default_rng(sample_seed)
            selected_ids = rng.choice(candidate_ids, size=max_proteins, replace=False)
        else:
            selected_ids = candidate_ids
        chosen_ids = set(str(value) for value in selected_ids)
    else:
        chosen_ids = target_protein_ids

    matched_records = []
    for record in stream_uniprot_fasta(fasta_path):
        if record.entry_id not in chosen_ids:
            continue
        if sample_seed is None and max_proteins is not None and len(matched_records) >= max_proteins:
            break
        matched_records.append(
            {
                "EntryID": record.entry_id,
                "entry_name": record.entry_name,
                "organism": record.organism,
                "sequence": record.sequence,
                "length": record.length,
                "terms": terms_grouped[record.entry_id],
            }
        )

    df = pd.DataFrame(
        matched_records,
        columns=["EntryID", "entry_name", "organism", "sequence", "length", "terms"],
    )
    return df
