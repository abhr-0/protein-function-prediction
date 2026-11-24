# ROADMAP.md — Milestones, Locked Decisions, Measured Facts

> Durable phase-level view. Last updated: 2026-11-24.
> Branch: `rebuild` (orphan, fresh history). `main` is the previous attempt and
> is reference only.
> `PLAN.md` is an optional local (untracked) commit-level checklist; `AGENTS.md`
> holds operating rules.

---

## Status at a glance

| Phase | Scope | Status |
|---|---|---|
| 0 | Foundation: docs, environment, raw data | **Complete (2026-11-24)** |
| 1 | Data adapters: human Swiss-Prot, CAFA6, evidence policy | Not started |
| 2 | Splits and vocabulary: MMseqs2 clusters, train-only vocab | Not started |
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
- [x] Counts verified: 20,429 proteins, 243,091 annotations; CAFA6 82,403
      train proteins / 537,027 annotations / 8,024 GT.

Acceptance met: `uv sync` works from a clean checkout, raw files match the
documented counts, and no code from the old branch is present. Phase 1 begins
with `feat: data adapters`.

---

## Phase 1 — Data adapters

- Human FASTA streaming + header parsing (accession, taxid, sequence).
- `DR GO` parsing from the derived human annotation file (and the REST txt
  stream); key on the primary accession; apply the evidence filter.
- CAFA6 adapters: `train_terms.tsv`, `train_taxonomy.tsv`, eval ground truth,
  `test_sequences.fasta`.
- Caches under `data/derived/cache/` (TSV/TSV.GZ; no parquet dependency).
- Tests: header parse, GO line parse, evidence policy, human filter, CAFA6
  counts.

Acceptance: human files yield 20,429 proteins and 243,091 annotations; the
`cafa` policy yields 216,272 annotations over 18,630 proteins and 16,442 terms;
CAFA6 train 82,403 proteins / 537,027 annotations; ground truth 8,024
(NK 749 / LK 485 / PK 6,951). No full flat file is read.

---

## Phase 2 — Splits and vocabulary

- MMseqs2 `easy-cluster` at 30% identity; split clusters 80/10/10.
- Train-only vocabulary; minimum training count policy; apply to val/test.
- Coverage report: terms kept, test annotations covered, proteins dropped with
  zero positives.
- Tests: no cluster spans splits; vocabulary is training-only; split is
  deterministic.

Acceptance: split JSON + vocabulary JSON + coverage report; no test protein
contributes to the vocabulary; each selected term has at least one train
positive.

---

## Phase 3 — ESM-2 embeddings

- `Embedder` protocol + ESM-2 650M implementation; mean-pool, truncate 1024.
- Cache key: model id + revision + FASTA SHA-256 + pooling + max length.
- Colab T4 notebook: mount Drive, upload bundle, fp16, sharded/resumable.
- Tests: cache round-trip with a tiny fake embedder; the key changes when the
  configuration changes.

Acceptance: embeddings for all 20,429 human proteins; deterministic rerun; cache
documented; Track B embedding plan sized (82,403 + 8,024 sequences).

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

- Train on CAFA6 t0 labels (82,403); predict the 8,024 ground-truth proteins.
- Prediction writer in CAFA format; `cafaeval` with the frozen graph + IA.
- PK uses `-known`; NK/LK as provided; per-aspect reporting.
- Tests: writer format; cafaeval wrapper on toy data.

Acceptance: reproducible cafaeval run for NK/LK/PK with the frozen
`go-basic.obo` and `IA.tsv`; outputs stored under `artifacts/runs/`; protocol
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
| Track A labels | UniProt `DR GO`, policy `cafa` (no IEA/ND/NR) |
| Track A split | MMseqs2 30% clusters, 80/10/10; no random comparison |
| Track A metric | Micro/macro P/R/F1; no IA/OBO in v1 |
| Track B universe | Full CAFA6: 82,403 train, 8,024 GT, 1,381 species |
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
| CAFA6 train | 82,403 proteins; 537,027 annotations; 1,381 species | `train_terms.tsv`, `train_taxonomy.tsv` |
| CAFA6 ground truth | 8,024 proteins: NK 749, LK 485, PK 6,951 | `May2026_groundtruth_*.tsv` |
| GT proteins in CAFA6 t0 train | 7,219 of 8,024 | ID intersection |
| GT proteins in Swiss-Prot 2026_04 | 8,023 of 8,024 | ID intersection |
| Human annotations (all evidence) | 243,091 over 17,317 terms; 19,375 proteins | REST txt stream |
| Human annotations (`cafa`) | 216,272 over 18,630 proteins; 16,442 terms | evidence filter |
| Human aspects (`cafa`) | F 46,519 / P 96,961 / C 72,792 | evidence filter |
| Human terms absent from frozen graph | 72 terms; 141 annotations (0.06%) | OBO vs annotations |
| Raw data on disk (gitignored) | 7.4 MiB UniProt + 295 MiB CAFA6 | `du -sh` |

---

## Why the frozen IA/graph

- `cafaeval` needs `-ia` for IA-weighted F-max; without it only unweighted F is
  produced.
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

- Track A is a contemporary-label, clustered holdout: not temporal, not
  CAFA-comparable.
- Clustering at 30% still leaves sub-threshold homology; numbers remain
  optimistic relative to a truly novel protein.
- Track B is the CAFA-comparable number, but PK scores only post-t0 terms and
  NK/LK are smaller.
- The frozen graph drops 0.06% of human annotations (no effect on Track A v1).
- ESM-2 is frozen; no fine-tuning in v1.
- Colab sessions are ephemeral; sharding is mandatory.
- No PPI, no external knowledge, no submission format in v1.
