"""Tests for the baseline runner's CLI defaults and warning logic.

The runner is a script, so it is loaded through importlib rather than imported
as a module. Loading has no side effects: work happens only inside ``main()``.
"""

import importlib.util
import sys
from argparse import Namespace
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_baseline_subset.py"


@pytest.fixture(scope="module")
def runner():
    """Load scripts/run_baseline_subset.py as a module."""
    spec = importlib.util.spec_from_file_location("run_baseline_subset", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_parse_args_defaults(runner, monkeypatch):
    """Default vocabulary policy must be count-based, not an arbitrary cap."""
    monkeypatch.setattr(sys, "argv", ["run_baseline_subset.py"])
    args = runner.parse_args()

    assert args.label_source == "uniprot"
    assert args.top_terms == 0, "a fixed top-50 cap silently dropped 60% of training data"
    assert args.min_term_count == 20
    assert args.min_terms == 1
    assert args.seed == 42
    assert args.no_sampling is False
    assert args.evidence_policy is None, "None distinguishes 'unset' from 'explicitly cafa'"


def test_environment_facts_are_not_command_line_flags(runner, monkeypatch):
    """Release and cache location are derived from provenance, never supplied.

    A user-supplied release could disagree with the bytes actually downloaded,
    which is exactly what the provenance manifest exists to prevent.
    """
    monkeypatch.setattr(sys, "argv", ["run_baseline_subset.py"])
    args = runner.parse_args()

    for removed in ("release", "cache_path", "sample_seed", "random_state"):
        assert not hasattr(args, removed), f"--{removed.replace('_', '-')} should not exist"


def test_detect_release_reads_manifest(runner, tmp_path, monkeypatch):
    """The release comes from the manifest, with a recorded source."""
    monkeypatch.setattr(
        runner, "get_data_manifest_path", lambda: tmp_path / "missing.json"
    )
    release, source = runner.detect_release()
    assert release == "unknown"
    assert source == "not recorded"


def _args(**overrides) -> Namespace:
    base = {
        "label_source": "uniprot",
        "n_proteins": 5000,
        "seed": 42,
        "no_sampling": False,
        "release": "2026_04",
    }
    base.update(overrides)
    return Namespace(**base)


def test_no_warnings_for_sane_configuration(runner):
    """The default configuration should be silent."""
    assert runner.build_warnings(_args(), explicit_evidence_policy=False) == []


def test_warns_when_evidence_policy_is_ignored(runner):
    """--evidence-policy has no meaning for CAFA labels."""
    warnings = runner.build_warnings(
        _args(label_source="cafa"), explicit_evidence_policy=True
    )
    assert any("evidence-policy" in warning for warning in warnings)


def test_warns_on_full_corpus_and_disabled_sampling(runner):
    """Both escape hatches are expensive or biased, so they must be announced."""
    warnings = runner.build_warnings(
        _args(n_proteins=0, no_sampling=True), explicit_evidence_policy=False
    )
    assert any("FULL Swiss-Prot corpus" in warning for warning in warnings)
    assert any("FASTA order" in warning for warning in warnings)


def test_warns_when_release_is_unknown(runner):
    """Unknown provenance must be visible, not silent."""
    warnings = runner.build_warnings(_args(release="unknown"), explicit_evidence_policy=False)
    assert any("release is unknown" in warning for warning in warnings)
