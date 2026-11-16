"""Data loading and preprocessing utilities for UniProt and CAFA."""

from src.data.annotations import (
    GoAnnotation,
    load_uniprot_go_annotations,
    parse_uniprot_go_line,
    resolve_excluded_evidence,
    stream_uniprot_go_annotations,
    summarize_annotations,
)
from src.data.uniprot import (
    UniProtRecord,
    create_annotated_dataset,
    load_annotations,
    load_information_accretion,
    parse_uniprot_header,
    stream_uniprot_fasta,
)

__all__ = [
    "UniProtRecord",
    "parse_uniprot_header",
    "stream_uniprot_fasta",
    "load_annotations",
    "load_information_accretion",
    "create_annotated_dataset",
    "GoAnnotation",
    "parse_uniprot_go_line",
    "stream_uniprot_go_annotations",
    "load_uniprot_go_annotations",
    "resolve_excluded_evidence",
    "summarize_annotations",
]
