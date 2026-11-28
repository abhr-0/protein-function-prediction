# Protein Function Prediction

Sequence-only prediction of Gene Ontology terms, rebuilt from scratch as two
experiments that share one codebase.

**Status: Phase 1 complete (2026-11-25). Phase 2 in progress: Task A complete
with the leakage-checked split; Task B (train-only vocabulary) in implementation.** The `pfpred`
package has adapters for both tracks, a `pfpred verify` command that checks the
25 raw-data counts, and `pfpred cluster` for Track A splits. `pfpred vocab` builds the per-aspect train-only GO vocabulary and coverage report (minimum five distinct training proteins per term).
See [ROADMAP.md](ROADMAP.md) for milestones. A local, untracked `PLAN.md` may
hold the current commit-level checklist.

## Two tracks

| | Track A — internal | Track B — CAFA6 |
|---|---|---|
| Proteins | human Swiss-Prot 2026_04 (20,429) | full CAFA6 (82,404 train / 8,023 GT) |
| Labels | UniProt `DR GO`, policy `curated` | CAFA6 `train_terms.tsv` (t0) |
| Split | MMseqs2 30% whole clusters; holdout-vs-train repair; ~83/8/8 | temporal: t0 → May 2026 |
| Metric | micro/macro P/R/F1 | cafaeval IA F-max, NK/LK/PK |
| CAFA-comparable | no | yes |

Aspects are trained separately, MF first. Features are frozen ESM-2 650M
embeddings; logistic regression and cosine kNN first, MLP later.

## Data

Track A uses two small files derived once from the UniProtKB REST stream
(`reviewed:true AND organism_id:9606`, release 2026_04):

```
data/raw/uniprot/human_proteins_2026_04.fasta.gz          # 20,429 proteins, 6.4 MB
data/raw/uniprot/human_go_annotations_2026_04.tsv.gz     # 243,091 rows: Entry, GO, Aspect, Evidence
```

Track B uses the full CAFA6 bundle:

```
data/raw/cafa6/train/{train_terms.tsv,train_sequences.fasta,train_taxonomy.tsv}
data/raw/cafa6/test/{testsuperset.fasta,testsuperset-taxon-list.tsv}
data/raw/cafa6/eval/{go-basic.obo,IA.tsv,test_sequences.fasta,
                     May2026_groundtruth_*.tsv,
                     May2026_groundtruth_terms_of_interest.txt,README.md}
data/raw/cafa6/sample_submission.tsv
```

`data/raw/` is read-only. `data/` and `artifacts/` are gitignored. There are no
download scripts; the UniProt query, derivation rules, hashes, and counts are in
[docs/data_sources.md](docs/data_sources.md).

## Environment

Nix + direnv provide Python, uv, and MMseqs2:

```bash
direnv allow        # or: nix develop
uv sync
uv run pytest
uv run pfpred verify
```

Python dependencies are managed by `uv`; see `pyproject.toml`. ESM-2 inference
is optional (`uv sync --extra embed`) and normally runs on Colab.

## Layout

```
configs/            # one YAML per track (later phases)
src/pfpred/         # adapters, policy, and CLI
tests/              # tiny fixtures, no network
notebooks/          # only the Colab embedding notebook
docs/               # data sources and methodology
data/               # gitignored: raw / derived
artifacts/          # gitignored: runs/<run_id>/
```

## Non-goals for v1

BP/CC, MLP, ProtT5, hand-crafted features, ESM-2 fine-tuning, IA-weighted F-max
on Track A, PPI, and download scripts are deliberately deferred.
