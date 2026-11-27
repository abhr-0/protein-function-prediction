# ROADMAP.md — Milestones, Locked Decisions, Measured Facts

> Durable phase-level view. Last updated: 2026-11-27.
> Branch: `rebuild` (orphan, fresh history). `main` is the previous attempt and
> is reference only.
> `PLAN.md` is an optional local (untracked) commit-level checklist; `AGENTS.md`
> holds operating rules.

---

## Status at a glance

| Phase | Scope | Status |
|---|---|---|
| 0 | Foundation: docs, environment, raw data | **Complete (2026-11-24)** |
| 1 | Data adapters: human Swiss-Prot, CAFA6, evidence policy | **Complete (2026-11-25)** |
| 2 | Splits and vocabulary: MMseqs2 clusters, train-only vocab | **In progress: Task A complete (2026-11-27); Task B vocabulary remains** |
| 3 | ESM-2 embeddings: protocol, cache, Colab T4 | Not started |
| 4 | Track A baseline: LR + kNN, MF, micro/macro F1 | Not started |
| 5 | Track B: full CAFA6, cafaeval, NK/LK/PK | Not started |
| 6 | Extensions: MLP, BP/CC, ProtT5, hand-crafted, fine-tuning | Not started |

---

## Phase 0 — Foundation (complete, 2026-11-24)

- [x] Orphan `rebuild` branch with fresh history.
- [x] `AGENTS.md`, `ROADMAP.md`, `README.md` (plus a local, untracked
      `PLAN.md`).
- [x] `pyproject.toml`, `flake.nix`, `.envrc`, `.gitignore`.
- [x] Place raw data in `data/raw/` (~302 MiB, gitignored and untracked):
      REST-derived UniProt files (7.4 MiB) + full CAFA6 bundle (295 MiB).
- [x] `direnv allow`/`nix develop` + `uv sync` succeeds; `numpy`, `sklearn`,
      `cafaeval` import.
- [x] Dev tooling: pytest 9.2.0, ruff 0.17.6.
- [x] Environment: Python 3.12.16, uv 0.11.21, numpy 2.5.5, scipy 1.18.1,
      scikit-learn 1.9.1, pandas 3.1.1, cafaeval 1.3.0; `uv.lock` tracked.
- [x] Raw counts verified: 20,429 proteins, 243,091 annotations; CAFA6 82,404
      train proteins / 537,027 annotations / 8,023 GT.

Acceptance met: `uv sync` works from a clean checkout, raw files match the
documented counts, and no code from the old branch is present. Phase 1 followed
with `build: make pfpred an installable package`.

---

## Phase 1 — Data adapters (complete, 2026-11-25)

- [x] `pfpred` packaged with hatchling; `pfpred` console script.
- [x] `records.py`: `Protein`, `Annotation`, `DataFormatError`, aspect validation.
- [x] `paths.py`: `Layout` with `PFP_ROOT` override, data paths, `sha256_file`.
- [x] `io.py`: gzip-aware text open, TSV rows, header validation.
- [x] `uniprot.py`: streaming human FASTA, annotation TSV, per-aspect label loading.
- [x] `labels.py`: `curated` evidence policy, grouping, summary statistics.
- [x] `cafa6.py`: CAFA6 train terms, headerless taxonomy, ground truth, sequences.
- [x] `cli.py`: `pfpred verify` checks every documented count and exits 1 on drift.
- [x] Tests: 81 tests on tiny fixtures; no network, no full-data reads.

Acceptance met: `pfpred verify` reports 20,429 proteins and 243,091 annotations;
the `curated` policy keeps 216,272 annotations over 18,630 proteins and 16,442
terms; CAFA6 train 82,404 proteins / 537,027 annotations; ground truth 8,023
(NK 749 / LK 485 / PK 6,951; 162 in both LK and PK). No full flat file is read.

---

## Phase 2 — Splits and vocabulary

### Task A — MMseqs2 clustering and protein splits

- [x] `pfpred cluster` command and `src/pfpred/splits.py` implementation.
- [x] MMseqs2 `easy-cluster`: 30% minimum identity, 80% bidirectional coverage
      (`-c 0.80 --cov-mode 0)).
- [x] Validate that every FASTA accession appears exactly once in cluster output.
- [x] Deterministic largest-cluster-first assignment; keep each cluster indivisible.
- [x] Record FASTA SHA-256, MMseqs2 version/parameters, seed, protein counts,
      cluster counts, and per-protein assignments in JSON.
- [x] Synthetic tests for TSV parsing, malformed/duplicate/missing members,
      deterministic assignment, and cluster integrity.
- [x] Run on the supplied Human Swiss-Prot FASTA; generated `data/derived/human_split.json`.
- [x] Inspect actual split counts: 20,429 proteins across 13,070 clusters; train 16,343 proteins / 10,461 clusters, validation 2,043 / 1,304, test 2,043 / 1,305.
- [x] Independently verify JSON integrity: metadata counts match assignments, all split labels are valid, and zero clusters span multiple splits.
- [x] Runtime validation on Windows: `pfpred verify`, `pytest`, and `ruff check .` all pass.

### Task B — Train-only vocabulary

- Train-only vocabulary; minimum training count policy; apply to val/test.
- Coverage report: terms kept, test annotations covered, proteins dropped with
  zero positives.
- Tests: vocabulary is training-only; selected terms have at least one train positive.

Acceptance for Phase 2: split JSON + vocabulary JSON + coverage report; no test
protein contributes to the vocabulary; each selected term has at least one train
positive. Phase 2 remains in progress until Task B is complete.

---

## Phase 3 — ESM-2 embeddings

- `Embedder` protocol + ESM-2 650M implementation; mean-pool, truncate 1024.
- Cache key: model id + revision + FASTA SHA-256 + pooling + max length.
- Colab T4 notebook: mount Drive, upload bundle, fp16, sharded/resumable.
- Tests: cache round-trip with a tiny fake embedder; the key changes when the
  configuration changes.

Acceptance: embeddings for all 20,429 human proteins; deterministic rerun; cache
documented; Track B embedding plan sized (82,404 + 8,023 sequences).

---

## Phase 4 — Track A baseline

- Logistic regression (binary relevance) and cosine kNN on ESM-2 features.
- MF first; validation-tuned threshold; test micro/macro P/R/F1.
- `artifacts/runs/<run_id>/` with config snapshot, metrics JSON, provenance.
- Tests: metric on toy predictions; no test leakage.

Acceptance: MF run completes end to end; metrics JSON records micro/macro P/R/F1,
threshold, coverage, and vocabulary size; Track A is labelled internal and
non-CAFA.

---

## Phase 5 — Track B (full CAFA6)

- Train on CAFA6 t0 labels (82,404); predict the 8,023 ground-truth proteins.
- Prediction writer in CAFA format; `cafaeval` with the frozen graph + IA.
- PK uses `-known`; NK/LK as provided; per-aspect reporting.
- Tests: writer format; cafaeval wrapper on toy data.

Acceptance: reproducible cafaeval run for NK/LK/PK with the frozen
`go-basic.obo` and `IA.tsv`, outputs stored under `artifacts/runs/`; protocol
flags recorded.

---

## Phase 6 — Extensions (after Phases 4–5)

- MLP (1–2 layers, BCE, early stopping).
- BP and CC aspects.
- ProtT5 and hand-crafted feature arms behind the same `Embedder` protocol.
- IA-weighted F-max for Track A (requires a graph/IA decision).
- Ancestor-propagated training labels (ablation).
- ESM-2 LoRA fine-tuning on Colab.
- Sequence-similarity kNN with MMseqs2 (homology baseline).
- PPI / external-knowledge track (explicitly out of scope until then).

---

## Locked decisions

| Decision | Choice |
|---|---|
| Branch | Orphan `rebuild`, fresh history |
| Environment | Nix flake + direnv; uv for Python dependencies |
| Data layout | `data/raw` (read-only inputs), `data/derived` (splits/vocab/embeddings), `artifacts/runs` |
| Track A universe | Human Swiss-Prot 2026_04, taxid 9606 |
| Track A labels | UniProt `DR GO`, evidence policy `curated` (no IEA/ND/NR) |
| Track A split | MMseqs2 `easy-cluster`, 30% identity, 80% bidirectional coverage (`-c 0.80 --cov-mode 0`), target 80/10/10 by protein count with whole clusters |
| Track A metric | Micro/macro P/R/F1; no IA/OBO in v1 |
| Track B universe | Full CAFA6: 82,404 train, 8,023 GT, 1,381 species |
| Track B labels | CAFA6 `train_terms.tsv` (t0) only |
| Track B metric | `cafaeval` IA F-max; NK/LK/PK; PK `-known` |
| GO graph/IA | CAFA6 frozen `go-basic.obo` + `IA.tsv` (Track B) |
| Propagation | Raw training labels; evaluate with `-prop fill` |
| Aspects | Per-aspect runs; MF first |
| Features | ESM-2 650M frozen, mean-pool, truncate 1024 |
| Models | Logistic regression + cosine kNN; MLP next |
| Embedding compute | Colab T4, fp16, sharded |
| Notebooks | Only `notebooks/colab_embed.ipynb` |
| Downloads | None; raw data copied manually |

---

## Measured facts

| Quantity | Value | How measured |
|---|---|---|
| Human Swiss-Prot 2026_04 | 20,429 proteins; 11,417,710 residues; mean 559; max 34,350 | FASTA scan |
| CAFA6 train | 82,404 proteins; 537,027 annotations; 1,381 species | `train_terms.tsv`, `train_taxonomy.tsv` |
| CAFA6 ground truth | 8,023 proteins: NK 749, LK 485, PK 6,951; 162 in both LK and PK | `May2026_groundtruth_*.tsv` |
| GT proteins in CAFA6 t0 train | 7,219 of 8,023 | ID intersection |
| GT proteins in Swiss-Prot 2026_04 | all 8,023; 2,232 are human | ID intersection |
| Human annotations (all evidence) | 243,091 over 17,317 terms; 19,375 proteins | REST txt stream |
| Human annotations (`curated`) | 216,272 over 18,630 proteins; 16,442 terms | evidence filter |
| Human aspects (`curated`) | F 46,519 / P 96,961 / C 72,792 | evidence filter |
| Human terms absent from frozen graph | 72 terms; 141 annotations (0.06%) | OBO vs annotations |
| Raw data on disk (gitignored) | 7.4 MiB UniProt + 295 MiB CAFA6 | `du -sh` |
| Verification | `uv run pfpred verify` checks 25 counts in ~5 s | exit 0 |
| Track A clustered split | 20,429 proteins in 13,070 clusters; train 16,343 / 10,461 clusters, validation 2,043 / 1,304, test 2,043 / 1,305 | `human_split.json`, seed 42 |
| Split integrity | 0 cross-split clusters; assignment and metadata counts match | Independent JSON integrity check |

---

## Why the frozen IA/graph

- `cafaeval` needs `-ia` for IA-weighted F-max; without it only unweighted F
  is produced.
- IA is tied to a specific OBO release. The only official pair available is
  CAFA6's 2025-06 `go-basic.obo` + `IA.tsv`.
- No official IA exists for "GO as of UniProt 2026_04".
- Track A v1 does not use IA or the graph.
- If Track A later gains IA F-max, using the frozen pair costs 72 of 17,317
  human terms and 141 of 243,091 annotations (0.06%); otherwise compute IA from
  the training labels (`IA(t) = -log2 P(t | parents)`, roughly 150 lines with
  `obonet`).

---

## Known limitations

- Track A uses contemporary labels and a clustered holdout; it is not a temporal
  or CAFA benchmark.
- Clustering at 30% identity and 80% bidirectional coverage still leaves
  sub-threshold homology; numbers remain optimistic relative to a truly novel protein.
- Track B is the CAFA-comparable number, but PK scores only post-t0 terms and
  NK/LK are smaller.
- The frozen graph drops 0.06% of human annotations (no effect on Track A v1).
- ESM-2 is frozen; no fine-tuning in v1.
- Colab sessions are ephemeral; sharding is mandatory.
