"""Tests for UniProt data streaming, header parsing, and annotation loading."""

import gzip
from pathlib import Path

import pandas as pd
import pytest

from src.data.uniprot import (
    UniProtRecord,
    build_annotated_dataset,
    load_annotations,
    load_information_accretion,
    normalize_aspect,
    parse_uniprot_header,
    stream_uniprot_fasta,
)
from src.utils.paths import get_ia_path, get_train_terms_path, get_uniprot_path


def test_parse_uniprot_header_full():
    """Verify regex parsing on standard full UniProt header."""
    header = (
        ">sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R "
        "OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-001R PE=4 SV=1"
    )
    parsed = parse_uniprot_header(header)

    assert parsed["entry_id"] == "Q6GZX4"
    assert parsed["entry_name"] == "001R_FRG3G"
    assert parsed["description"] == "Putative transcription factor 001R"
    assert parsed["organism"] == "Frog virus 3 (isolate Goorha)"
    assert parsed["taxon_id"] == 654924
    assert parsed["gene_name"] == "FV3-001R"
    assert parsed["protein_evidence"] == 4
    assert parsed["sequence_version"] == 1


def test_parse_uniprot_header_without_gene():
    """Verify parsing when optional fields like GN are absent."""
    header = ">sp|Q6GZX3|002L_FRG3G Uncharacterized protein 002L OS=Frog virus 3 OX=654924 PE=4 SV=1"
    parsed = parse_uniprot_header(header)

    assert parsed["entry_id"] == "Q6GZX3"
    assert parsed["entry_name"] == "002L_FRG3G"
    assert parsed["gene_name"] is None
    assert parsed["taxon_id"] == 654924


def test_normalize_aspect():
    """Verify normalization of various aspect string representations."""
    assert normalize_aspect("F") == "F"
    assert normalize_aspect("mf") == "F"
    assert normalize_aspect("molecular_function") == "F"
    assert normalize_aspect("P") == "P"
    assert normalize_aspect("C") == "C"

    with pytest.raises(ValueError):
        normalize_aspect("INVALID")


def test_stream_uniprot_fasta_synthetic(tmp_path: Path):
    """Test streaming from a temporary gzip-compressed FASTA file."""
    fasta_gz = tmp_path / "sample.fasta.gz"
    content = (
        ">sp|P12345|TEST1_HUMAN Test protein 1 OS=Homo sapiens OX=9606 GN=TST1 PE=1 SV=1\n"
        "MKTIIALSYIFCLVFADYKDDDDK\n"
        ">sp|P67890|TEST2_HUMAN Test protein 2 OS=Homo sapiens OX=9606 GN=TST2 PE=1 SV=1\n"
        "MALLALALALALALALALAL\n"
    )

    with gzip.open(fasta_gz, "wt", encoding="utf-8") as f:
        f.write(content)

    records = list(stream_uniprot_fasta(fasta_gz))
    assert len(records) == 2
    assert isinstance(records[0], UniProtRecord)
    assert records[0].entry_id == "P12345"
    assert records[0].sequence == "MKTIIALSYIFCLVFADYKDDDDK"
    assert records[0].length == 24
    assert records[1].entry_id == "P67890"

    # Test max_records limit
    records_limited = list(stream_uniprot_fasta(fasta_gz, max_records=1))
    assert len(records_limited) == 1


def test_load_annotations_synthetic(tmp_path: Path):
    """Test annotation loading and filtering with synthetic TSV."""
    terms_file = tmp_path / "terms.tsv"
    data = (
        "EntryID\tterm\taspect\n"
        "P1\tGO:0001\tF\n"
        "P1\tGO:0002\tF\n"
        "P2\tGO:0001\tF\n"
        "P2\tGO:0003\tP\n"
        "P3\tGO:0004\tC\n"
    )
    terms_file.write_text(data)

    # Filter for MF (F)
    df = load_annotations(terms_file, aspect="F")
    assert len(df) == 3
    assert set(df["term"]) == {"GO:0001", "GO:0002"}

    # No aspect filter keeps every row
    df_all = load_annotations(terms_file)
    assert len(df_all) == 5
    assert set(df_all["aspect"]) == {"F", "P", "C"}


def test_production_path_load_then_build(tmp_path: Path):
    """CAFA labels flow through load_annotations -> build_annotated_dataset."""
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(
        ">sp|P1|PROT1_HUMAN Protein 1 OS=Homo sapiens OX=9606 PE=1 SV=1\n"
        "ACDEFGHIKLMNPQRSTVWY\n"
        ">sp|P2|PROT2_HUMAN Protein 2 OS=Homo sapiens OX=9606 PE=1 SV=1\n"
        "YYYYYYYYYY\n"
        ">sp|P3|PROT3_HUMAN Unannotated Protein OS=Homo sapiens OX=9606 PE=1 SV=1\n"
        "AAAAAAAAAA\n"
    )

    terms_file = tmp_path / "terms.tsv"
    terms_file.write_text(
        "EntryID\tterm\taspect\n"
        "P1\tGO:0001\tF\n"
        "P1\tGO:0002\tF\n"
        "P2\tGO:0001\tF\n"
        "P2\tGO:0003\tP\n"
        "P2\tGO:0004\tC\n"
    )

    # aspect=None keeps all aspects
    annotations = load_annotations(terms_file)
    df = build_annotated_dataset(fasta_file, annotations)
    assert set(df["EntryID"]) == {"P1", "P2"}  # P3 has no annotations
    p2_terms = df[df["EntryID"] == "P2"].iloc[0]["terms"]
    assert p2_terms == ["GO:0001", "GO:0003", "GO:0004"]

    # aspect filter drops the other ontologies
    annotations_f = load_annotations(terms_file, aspect="F")
    df_f = build_annotated_dataset(fasta_file, annotations_f, aspect="F")
    p2_terms_f = df_f[df_f["EntryID"] == "P2"].iloc[0]["terms"]
    assert p2_terms_f == ["GO:0001"]


def test_real_dataset_smoke():
    """Smoke test: verify reading real workspace files (small slice)."""
    uniprot_path = get_uniprot_path()
    terms_path = get_train_terms_path()
    ia_path = get_ia_path()

    if not (uniprot_path.exists() and terms_path.exists() and ia_path.exists()):
        pytest.skip("Workspace raw data files not found.")

    # Stream 5 real records
    records = list(stream_uniprot_fasta(uniprot_path, max_records=5))
    assert len(records) == 5
    for r in records:
        assert len(r.entry_id) > 0
        assert len(r.sequence) > 0
        assert r.length == len(r.sequence)

    # Load IA weights
    ia_dict = load_information_accretion(ia_path)
    assert len(ia_dict) > 1000
    assert all(isinstance(v, float) for v in list(ia_dict.values())[:10])

    # Create tiny dataset of 10 proteins through the production path
    annotations = load_annotations(terms_path, aspect="F")
    df = build_annotated_dataset(
        uniprot_path,
        annotations,
        aspect="F",
        max_proteins=10,
    )
    assert len(df) == 10
    assert "sequence" in df.columns
    assert "terms" in df.columns


def _synthetic_annotations(entries: list[str]) -> pd.DataFrame:
    """Build a tidy annotation frame for the generic builder tests."""
    rows = []
    for index, entry_id in enumerate(entries):
        rows.append({"EntryID": entry_id, "term": "GO:0001", "aspect": "F"})
        if index % 2 == 0:
            rows.append({"EntryID": entry_id, "term": "GO:0002", "aspect": "F"})
    return pd.DataFrame(rows)


def test_build_annotated_dataset_generic(tmp_path: Path):
    """The source-agnostic builder should align arbitrary tidy annotations."""
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(
        ">sp|P1|PROT1_HUMAN Protein 1 OS=Homo sapiens OX=9606 PE=1 SV=1\n"
        "ACDEFGHIKLMNPQRSTVWY\n"
        ">sp|P2|PROT2_HUMAN Protein 2 OS=Homo sapiens OX=9606 PE=1 SV=1\n"
        "YYYYYYYYYY\n"
        ">sp|P3|PROT3_HUMAN Protein 3 OS=Homo sapiens OX=9606 PE=1 SV=1\n"
        "AAAAAAAAAA\n"
    )
    annotations = pd.DataFrame(
        [
            {"EntryID": "P1", "term": "GO:0001", "aspect": "F"},
            {"EntryID": "P1", "term": "GO:0002", "aspect": "F"},
            {"EntryID": "P2", "term": "GO:0003", "aspect": "F"},
            {"EntryID": "P9", "term": "GO:0001", "aspect": "F"},
        ]
    )

    df = build_annotated_dataset(
        fasta_file,
        annotations,
        target_terms=["GO:0001", "GO:0002"],
    )

    assert set(df["EntryID"]) == {"P1"}  # P2 only has GO:0003; P9 is absent
    assert df.iloc[0]["terms"] == ["GO:0001", "GO:0002"]


def test_build_annotated_dataset_aspect_filter(tmp_path: Path):
    """The aspect filter should restrict annotations to one ontology."""
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(
        ">sp|P1|PROT1_HUMAN Protein 1 OS=Homo sapiens OX=9606 PE=1 SV=1\n"
        "ACDEFGHIKLMNPQRSTVWY\n"
        ">sp|P2|PROT2_HUMAN Protein 2 OS=Homo sapiens OX=9606 PE=1 SV=1\n"
        "YYYYYYYYYY\n"
    )
    annotations = pd.DataFrame(
        [
            {"EntryID": "P1", "term": "GO:0001", "aspect": "F"},
            {"EntryID": "P1", "term": "GO:0002", "aspect": "F"},
            {"EntryID": "P2", "term": "GO:0003", "aspect": "F"},
            {"EntryID": "P2", "term": "GO:0004", "aspect": "C"},
        ]
    )

    df = build_annotated_dataset(fasta_file, annotations, aspect="F")

    assert set(df["EntryID"]) == {"P1", "P2"}
    assert set(df[df["EntryID"] == "P2"].iloc[0]["terms"]) == {"GO:0003"}


def test_build_annotated_dataset_seeded_sampling_is_deterministic(tmp_path: Path):
    """Random protein sampling should be reproducible and bounded."""
    fasta_file = tmp_path / "test.fasta"
    lines = []
    for index in range(10):
        accession = f"P{index:05d}"
        lines.append(
            f">sp|{accession}|PROT{index}_HUMAN Protein {index} "
            "OS=Homo sapiens OX=9606 PE=1 SV=1\n"
            "ACDEFGHIKLMNPQRSTVWY\n"
        )
    fasta_file.write_text("".join(lines))

    annotations = _synthetic_annotations([f"P{index:05d}" for index in range(10)])

    first = build_annotated_dataset(
        fasta_file,
        annotations,
        max_proteins=5,
        sample_seed=42,
    )
    second = build_annotated_dataset(
        fasta_file,
        annotations,
        max_proteins=5,
        sample_seed=42,
    )

    assert len(first) == 5
    assert set(first["EntryID"]) == set(second["EntryID"])
    assert set(first["EntryID"]).issubset(set(annotations["EntryID"]))


def test_build_annotated_dataset_empty_annotations(tmp_path: Path):
    """Empty annotations should yield an empty frame with the documented schema."""
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(
        ">sp|P1|PROT1_HUMAN Protein 1 OS=Homo sapiens OX=9606 PE=1 SV=1\n"
        "ACDEFGHIKLMNPQRSTVWY\n"
    )
    annotations = pd.DataFrame(columns=["EntryID", "term", "aspect"])

    df = build_annotated_dataset(fasta_file, annotations)

    assert len(df) == 0
    assert list(df.columns) == [
        "EntryID",
        "entry_name",
        "organism",
        "sequence",
        "length",
        "terms",
    ]
