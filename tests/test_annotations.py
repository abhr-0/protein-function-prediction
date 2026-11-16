"""Tests for UniProt flat-file GO annotation parsing and caching."""

import gzip
from pathlib import Path

import pandas as pd
import pytest

from src.data.annotations import (
    EVIDENCE_POLICIES,
    load_uniprot_go_annotations,
    parse_uniprot_go_line,
    resolve_excluded_evidence,
    stream_uniprot_go_annotations,
    summarize_annotations,
)

SYNTHETIC_DAT = """\
ID   TEST1_HUMAN             Reviewed;         100 AA.
AC   P12345; Q99999;
DR   GO; GO:0004842; F:ubiquitin-protein transferase activity; IDA:UniProtKB.
DR   GO; GO:0005515; F:protein binding; IEA:UniProtKB-KW.
DR   GO; GO:0005737; C:cytoplasm; IDA:UniProtKB.
DR   GO; GO:0008150; P:biological_process; ND:UniProtKB.
//
ID   TEST2_HUMAN             Reviewed;         100 AA.
AC   P67890;
DR   GO; GO:0005524; F:ATP binding; IEA:UniProtKB.
DR   GO; GO:0005634; C:nucleus; IDA:UniProtKB.
DR   GO; GO:0006468; P:protein phosphorylation; IMP:UniProtKB.
//
ID   TEST3_HUMAN             Reviewed;         100 AA.
AC   Q11111;
//
"""


@pytest.fixture()
def synthetic_dat(tmp_path: Path) -> Path:
    """Write a small gzip-compressed UniProt flat file."""
    path = tmp_path / "synthetic.dat.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(SYNTHETIC_DAT)
    return path


def test_parse_uniprot_go_line():
    """Fields should be extracted without relying on fixed column widths."""
    line = (
        "DR   GO; GO:0004842; F:ubiquitin-protein transferase activity; "
        "IDA:UniProtKB.\n"
    )
    annotation = parse_uniprot_go_line(line)

    assert annotation is not None
    assert annotation.term == "GO:0004842"
    assert annotation.aspect == "F"
    assert annotation.evidence == "IDA"


def test_parse_uniprot_go_line_rejects_non_go_and_malformed():
    """Non-GO lines and junk should return None rather than raise."""
    assert parse_uniprot_go_line("DR   GO; not-a-go-line\n") is None
    assert parse_uniprot_go_line("DR   GO; GO:0001; X:bad aspect; IDA:Test.\n") is None
    assert parse_uniprot_go_line("DR   EM; 1; 2; 3.\n") is None


def test_resolve_excluded_evidence():
    """Named policies and explicit collections should both resolve."""
    assert resolve_excluded_evidence("cafa") == EVIDENCE_POLICIES["cafa"]
    assert resolve_excluded_evidence("FULL") == frozenset({"ND", "NR"})
    assert resolve_excluded_evidence(["IEA", "nd"]) == frozenset({"IEA", "ND"})
    with pytest.raises(ValueError):
        resolve_excluded_evidence("not-a-policy")


def test_stream_uses_primary_accession_and_filters(synthetic_dat: Path):
    """The first AC value should be the key; IEA/ND are excluded by default."""
    records = list(stream_uniprot_go_annotations(synthetic_dat))

    assert {r.entry_id for r in records} == {"P12345", "P67890"}
    assert all(r.entry_id == "P12345" or r.entry_id == "P67890" for r in records)
    assert ("P12345", "GO:0004842", "F", "IDA") in {
        (r.entry_id, r.term, r.aspect, r.evidence) for r in records
    }
    assert "GO:0005515" not in {r.term for r in records}  # IEA
    assert "GO:0008150" not in {r.term for r in records}  # ND
    assert "GO:0005524" not in {r.term for r in records}  # IEA


def test_stream_full_policy_includes_iea(synthetic_dat: Path):
    """The full policy keeps IEA but still drops ND/NR."""
    records = list(
        stream_uniprot_go_annotations(
            synthetic_dat,
            excluded_evidence=EVIDENCE_POLICIES["full"],
        )
    )
    terms = {r.term for r in records}

    assert "GO:0005515" in terms
    assert "GO:0005524" in terms
    assert "GO:0008150" not in terms


def test_stream_aspect_filter(synthetic_dat: Path):
    """Aspect filtering should happen during streaming."""
    records = list(stream_uniprot_go_annotations(synthetic_dat, aspects=["C"]))
    assert {r.aspect for r in records} == {"C"}
    assert {r.term for r in records} == {"GO:0005737", "GO:0005634"}


def test_load_returns_tidy_frame(synthetic_dat: Path):
    """The loader should return the documented schema with categoricals."""
    frame = load_uniprot_go_annotations(synthetic_dat, aspects=["F"])

    assert list(frame.columns) == ["EntryID", "term", "aspect", "evidence"]
    assert len(frame) == 1
    assert frame.iloc[0]["term"] == "GO:0004842"
    assert isinstance(frame["EntryID"].dtype, pd.CategoricalDtype)


def test_load_writes_and_reuses_cache(synthetic_dat: Path, tmp_path: Path):
    """The cache should be written once and reused on the next call."""
    cache = tmp_path / "cache" / "annotations.tsv.gz"

    first = load_uniprot_go_annotations(synthetic_dat, cache_path=cache)
    assert cache.is_file()
    assert len(first) == 4  # P12345 (2) + P67890 (2), IEA/ND excluded

    second = load_uniprot_go_annotations(synthetic_dat, cache_path=cache)
    pd.testing.assert_frame_equal(first, second)


def test_load_force_rebuilds_cache(synthetic_dat: Path, tmp_path: Path):
    """force=True should re-parse even if a cache exists."""
    cache = tmp_path / "annotations.tsv.gz"
    load_uniprot_go_annotations(synthetic_dat, cache_path=cache)

    cache.write_bytes(b"corrupt")
    frame = load_uniprot_go_annotations(synthetic_dat, cache_path=cache, force=True)

    assert len(frame) == 4


def test_load_missing_file_raises(tmp_path: Path):
    """Missing flat files should raise a clear FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_uniprot_go_annotations(tmp_path / "absent.dat.gz")


def test_summarize_annotations(synthetic_dat: Path):
    """Summary counts should match the synthetic flat file."""
    frame = load_uniprot_go_annotations(
        synthetic_dat,
        evidence_policy=EVIDENCE_POLICIES["full"],
    )
    summary = summarize_annotations(frame)

    assert summary["annotations"] == 6
    assert summary["proteins"] == 2
    assert summary["aspects"]["F"]["proteins"] == 2
    assert summary["aspects"]["C"]["annotations"] == 2
    assert summary["aspects"]["P"]["annotations"] == 1
    assert summary["evidence"]["IEA"] == 2
    assert summary["evidence"].get("ND", 0) == 0


def test_summarize_empty_annotations():
    """An empty frame should summarize without raising."""
    frame = pd.DataFrame(columns=["EntryID", "term", "aspect", "evidence"])
    summary = summarize_annotations(frame)
    assert summary == {
        "annotations": 0,
        "proteins": 0,
        "terms": 0,
        "aspects": {},
        "evidence": {},
    }



MULTI_AC_DAT = """\
ID   PEX5R_HUMAN             Reviewed;         626 AA.
AC   Q8IYB4; B7Z2A5; B7Z305; B7Z318; B7Z5Z5; B7Z8P2; E7EUV8; E7EUZ0; E9PEC1;
AC   E9PH97; Q9NQD1; Q9P2U3; Q9P2U4;
DR   GO; GO:0005515; F:protein binding; IPI:UniProtKB.
//
"""


def test_stream_uses_primary_accession_across_wrapped_ac_lines(tmp_path: Path):
    """Only the first value of the first AC line is the primary accession.

    Regression test: the AC field wraps for merged entries, and re-keying on
    every AC line attributed annotations to a secondary accession (e.g.
    E9PH97 instead of Q8IYB4), which is absent from the FASTA and therefore
    silently dropped during alignment.
    """
    path = tmp_path / "multi_ac.dat.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(MULTI_AC_DAT)

    records = list(stream_uniprot_go_annotations(path))
    frame = load_uniprot_go_annotations(path)

    assert [r.term for r in records] == ["GO:0005515"]
    assert {r.entry_id for r in records} == {"Q8IYB4"}
    assert "E9PH97" not in set(frame["EntryID"].astype(str))
