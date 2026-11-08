"""UniProtKB/Swiss-Prot data streaming, header parsing, and annotation joining."""

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

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
    min_term_count: int | None = None,
    top_k_terms: int | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """Load and optionally filter GO annotations from train_terms.tsv.

    Args:
        terms_path: Path to train_terms.tsv.
        aspect: Optional filter ('F', 'P', or 'C').
        min_term_count: Minimum frequency for a GO term to be included.
        top_k_terms: Top K most frequent GO terms to keep.

    Returns:
        tuple of (filtered_terms_dataframe, list_of_selected_terms)
    """
    terms_file = Path(terms_path)
    if not terms_file.exists():
        raise FileNotFoundError(f"Annotations file not found: {terms_file}")

    df = pd.read_csv(terms_file, sep="\t")

    if aspect is not None:
        normalized = normalize_aspect(aspect)
        df = df[df["aspect"] == normalized]

    term_counts = df["term"].value_counts()

    if min_term_count is not None:
        valid_terms = term_counts[term_counts >= min_term_count].index
        df = df[df["term"].isin(valid_terms)]
        term_counts = term_counts[term_counts.index.isin(valid_terms)]

    if top_k_terms is not None:
        top_terms = term_counts.head(top_k_terms).index
        df = df[df["term"].isin(top_terms)]
        selected_terms = list(top_terms)
    else:
        selected_terms = list(term_counts.index)

    return df, selected_terms


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


def create_annotated_dataset(
    fasta_path: str | Path,
    terms_path: str | Path,
    aspect: str | None = None,
    max_proteins: int | None = 5000,
    min_term_count: int | None = None,
    top_k_terms: int | None = 300,
) -> tuple[pd.DataFrame, list[str]]:
    """Build an aligned DataFrame of sequences and multi-label GO annotations.

    Streams UniProt FASTA records and aligns them with their associated Gene
    Ontology terms from a tabular annotations file (e.g. CAFA train_terms.tsv).

    Args:
        fasta_path: Path to plain or gzip-compressed FASTA file.
        terms_path: Path to annotations TSV file with columns (EntryID, term, aspect).
        aspect: Optional filter for a specific Gene Ontology sub-ontology:
            - 'F' or 'MF': Molecular Function (biochemical activity)
            - 'P' or 'BP': Biological Process (cellular pathway / program)
            - 'C' or 'CC': Cellular Component (subcellular localization)
            - None (default): Retains all aspects jointly. When using a limited
              `top_k_terms`, frequent Cellular Component terms may dominate the top
              ranks; callers seeking balanced representation across all categories
              can either set `top_k_terms=None` (filtering by `min_term_count`
              instead) or call this function separately for each aspect.
        max_proteins: Maximum number of annotated proteins to yield (None for all).
        min_term_count: Minimum occurrence frequency for a GO term to be included.
            Useful for pruning rare terms without aspect bias.
        top_k_terms: Top K most frequent GO terms to retain (None for all).

    Returns:
        tuple of (dataset_df, target_terms_list):
            - dataset_df: DataFrame with columns ['EntryID', 'entry_name',
              'organism', 'sequence', 'length', 'terms']
            - target_terms_list: Ordered list of unique target GO term IDs.
    """
    annotations_df, target_terms = load_annotations(
        terms_path=terms_path,
        aspect=aspect,
        min_term_count=min_term_count,
        top_k_terms=top_k_terms,
    )

    # Group terms per protein
    terms_grouped = (
        annotations_df.groupby("EntryID")["term"]
        .apply(lambda s: sorted(list(set(s))))
        .to_dict()
    )
    target_protein_ids = set(terms_grouped.keys())

    matched_records = []
    for record in stream_uniprot_fasta(fasta_path):
        if record.entry_id in target_protein_ids:
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
            if max_proteins is not None and len(matched_records) >= max_proteins:
                break

    df = pd.DataFrame(matched_records)
    return df, target_terms
