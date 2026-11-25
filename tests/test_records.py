from __future__ import annotations

import pytest

from pfpred.records import Annotation, DataFormatError, Protein, validate_aspect


@pytest.mark.parametrize("aspect", ["F", "P", "C"])
def test_validate_aspect_accepts_known_letters(aspect: str) -> None:
    assert validate_aspect(aspect, source="test") == aspect


@pytest.mark.parametrize("aspect", ["f", "MF", "", "X"])
def test_validate_aspect_rejects_unknown_letters(aspect: str) -> None:
    with pytest.raises(DataFormatError, match="unknown aspect"):
        validate_aspect(aspect, source="test")


def test_annotation_evidence_defaults_to_empty() -> None:
    annotation = Annotation(entry="P12345", term="GO:0003674", aspect="F")
    assert annotation.evidence == ""


def test_records_are_frozen() -> None:
    protein = Protein(accession="P12345", sequence="MKV")
    with pytest.raises(AttributeError):
        protein.sequence = "M"  # type: ignore[misc]
