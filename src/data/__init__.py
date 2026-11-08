"""Data loading and preprocessing utilities for UniProt and CAFA."""

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
]
