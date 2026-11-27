# AGENTS.md — Operating Manual for the Rebuild

> **Project:** sequence-only protein function prediction.
> **Branch:** `rebuild` (orphan branch; `main` is the previous attempt and is
> reference only).
> **Status:** Phase 1 complete (2026-11-25). Phase 2 is in progress: Task A
> clustering/splits completed and validated on real data (2026-11-27); Task B
> train-only vocabulary remains.
> **Companions:** `ROADMAP.md` (milestones, locked decisions, measured facts).
> `PLAN.md` is an optional local working plan, untracked and possibly stale.

---

## 1. Objective and scope

Given an amino-acid sequence, predict Gene Ontology (GO) terms for Molecular
Function (`F`), Biological Process (`P`), and Cellular Component (`C`). v1
focuses on MF.

The rebuild runs two experiments on one codebase:

| | Track A — internal benchmark | Track B — CAFA6 benchmark |
|---|---|---|
| Protein universe | human Swiss-Prot, taxid 9606, 20,429 proteins | full CAFA6: 82,404 train proteins, 1,381 species |
| Labels | UniProt `DR GO` (human REST txt, release 2026_04), evidence policy `curated` | CAFA6 `train/train_terms.tsv` (t0 snapshot) |
| Split | MMseqs2 30% identity + 80% bidirectional coverage, 80/10/10 by whole cluster | temporal: t0 labels → May 2026 ground truth |
| Evaluation | micro/macro precision, recall, F1 at a validation-tuned threshold | `cafaeval` IA-weighted F-max, NK/LK/PK |
| CAFA-comparable | no | yes (with documented caveats) |
| Aspect model | one run per aspect, MF first | one run per aspect, MF first |

**Do not merge the tracks.** They answer different questions with different label
vintages and metrics. A Track A F1 is never a CAFA result, and a Track B number
is never a Track A holdout result.

---

## 2. Measured facts (do not re-derive without cause)

| Quantity | Value |
|---|---|
| Human Swiss-Prot 2026_04 | 20,429 proteins; 11,417,710 residues; mean 559 aa; max 34,350 |
| CAFA6 train | 82,404 proteins; 537,027 annotations; 1,381 species |
| CAFA6 ground truth | 8,023 proteins: NK 749, LK 485, PK 6,951 (162 in both LK and PK) |
| GT proteins present in CAFA6 t0 train | 7,219 of 8,023 |
| GT proteins present in Swiss-Prot 2026_04 | all 8,023; 2,232 are human |
| Human annotations (all evidence) | 243,091 over 17,317 terms; 19,375 proteins |
| Human annotations under `curated` policy | 216,272 over 18,630 proteins; 16,442 terms |
| Human aspect counts (`curated`) | F 46,519 / P 96,961 / C 72,792 |
| Human terms absent from frozen 2025-06 graph | 72 terms; 141 annotations (0.06%) |
| Embedding compute | Colab T4, fp16, sharded |

---

## 3. Hard rules

### 3.1 Leakage

- Split proteins **before** vocabulary selection, feature fitting, or threshold
  selection.
- Build the term vocabulary on the **training split only**, then apply it to
  validation and test. Report test coverage and zero-positive drops.
- Track B labels must come from CAFA6 `train_terms.tsv` (t0). **Never** train
  Track B on UniProt 2026_04 labels: all 8,023 GT proteins are in the full
  Swiss-Prot 2026_04 flat file and 2,232 are human, so the human files can carry
  post-t0 annotations (7,219 GT proteins already had t0 labels).
- Track A and Track B each have their own vocabulary; never share one across
  tracks.
- Per-aspect protein sets: use proteins with at least one label in that aspect.
  Do not intersect MF/BP/CC protein sets.
- Test data never influences `fit`, vocabulary, thresholds, or hyperparameter
  choice.
- `Embedder.fit` (SVD, scalers) runs on training sequences only; frozen
  embedders implement it as a no-op.

### 3.2 Labels

- Evidence policy `curated`: exclude `IEA`, `ND`, `NR`; keep the rest. The name
  avoids colliding with Track B's CAFA6 benchmark.
- Train on raw annotations. Do not propagate to GO ancestors during training.
- Track A v1 reads no ontology: no OBO, no IA, no propagation. Its metric is F1.
- Track B propagates only at evaluation, via `cafaeval -prop fill`.

### 3.3 Evaluation

- Track A: micro/macro precision, recall, F1; one threshold per aspect tuned on
  validation, applied unchanged to test.
- Track B: `cafaeval` with CAFA6's frozen `go-basic.obo` + `IA.tsv`; official
  protocol: PK uses `-known`, plus
  `-no_orphans -norm cafa -prop fill -toi ... -max_terms 500 -th_step 0.001`.
- Report label-space coverage for every run.
- Label every number with its track, aspect, split, and metric. Never present one
  track's number as the other's.

### 3.4 Features

- v1 features are frozen ESM-2 650M (`facebook/esm2_t33_650M_UR50D`),
  mean-pooled over residues, truncated at 1024 aa.
- Every representation implements the `Embedder` protocol: `name`, `dim`,
  `fit(sequences)`, `transform(sequences) -> ndarray`, `cache_key()`.
- The cache key includes model id, revision, FASTA SHA-256, pooling, and max
  length. Never reuse a cache with a different key.
- Store pooled vectors only (~105 MB Track A, ~420 MB Track B). Per-residue
  storage is out of scope.
- No fine-tuning in v1.

### 3.5 Data

- `data/raw/` is read-only. Reusable derived data (splits, vocabulary,
  embeddings) goes to `data/derived/`; run outputs go to
  `artifacts/runs/<run_id>/`.
- There are **no download scripts**. UniProt data is derived once from the REST
  stream and documented in `docs/data_sources.md`; the CAFA6 bundle is copied
  from the competition/evaluation release.
- Track A inputs: `data/raw/uniprot/human_proteins_2026_04.fasta.gz` and
  `data/raw/uniprot/human_go_annotations_2026_04.tsv.gz`
  (`Entry`, `GO`, `Aspect`, `Evidence`).
- Track B inputs: the full CAFA6 train/test/eval bundle under `data/raw/cafa6/`.
- v1 does **not** need the full `uniprot_sprot.dat.gz`, the full Swiss-Prot
  FASTA, `data/go/go-basic.obo`, or a UniProt web TSV.

---

## 4. Environment and commands

Nix flake + direnv provide Python 3.12, uv, and MMseqs2; uv owns Python
dependencies.

```bash
direnv allow                 # or: nix develop
uv sync                      # core + dev dependencies
uv sync --extra embed        # torch + transformers (usually only on Colab)
uv run pytest
uv run ruff check .
uv run pfpred verify         # raw files vs documented counts
uv run pfpred cluster        # MMseqs2 clusters + Track A splits
```

Do not add dependencies without a reason recorded in `ROADMAP.md` (and, locally,
`PLAN.md`).
Keep the core set small: numpy, scipy, scikit-learn, pandas, pyyaml, cafaeval.

---

## 5. Planned layout

```
configs/             # one YAML per track
src/pfpred/
  records.py         # Protein, Annotation, formatting errors
  paths.py           # repository layout and data paths
  io.py              # gzip-aware text and TSV helpers
  uniprot.py         # human FASTA + annotation adapters
  cafa6.py           # CAFA6 train/GT adapters
  labels.py          # evidence policy, grouping, summary stats
  splits.py          # MMseqs2 clusters -> split
  features/
    base.py          # Embedder protocol + registry
    esm2.py          # ESM-2 extraction and cache
  models.py          # logistic regression, kNN
  evaluate.py        # Track A F1; Track B cafaeval wrapper
  cli.py             # verify | cluster | embed | train | evaluate
tests/               # tiny fixtures, no network
notebooks/           # only notebooks/colab_embed.ipynb
docs/                # data_sources.md, methodology.md
data/                # gitignored: raw/ (inputs), derived/ (splits/vocab/embeddings)
artifacts/           # gitignored: runs/<run_id>/
```

---

## 6. Working agreement

- `ROADMAP.md` is authoritative for phase status and locked decisions.
  `PLAN.md`, if present locally, is the current commit-level checklist; it is
  untracked and may be stale. If they disagree, follow `ROADMAP.md` and say so.
- Use **Conventional Commits** for every commit: `type(scope): summary`, with
  `feat`, `fix`, `docs`, `chore`, `test`, `refactor`, `perf`, or `build`
  as the type. Examples: `feat: data adapters`, `fix: key annotations on primary
  accession`. Do not use free-form commit messages.
- Work in small commits that match the current plan (`PLAN.md` if present
  locally). Meet a phase's acceptance criteria before starting the next.
- Every module gets a test against tiny synthetic fixtures. No network, no
  700 MB files, no real ESM-2 in tests.
- Use the standard library for FASTA and `.dat.gz` parsing; stream the flat file;
  never load it whole.
- Keep matrices sparse where the representation is sparse; dense is fine for
  20k × 1280 embeddings.
- No notebooks except the Colab embedding notebook. Plots and EDA, if any, live
  in `docs/` or `artifacts/`.
- Do not reintroduce the old `src/`, `scripts/download_data.py`, or the old
  CAFA-only pipeline.
- Adding BP/CC, MLP, ProtT5, hand-crafted features, fine-tuning, or PPI requires
  a `ROADMAP.md` update first.

---

## 7. Current status and next step

Phase 1 is complete: the `pfpred` package ships adapters for the human UniProt
files and the CAFA6 bundle, the `curated` evidence policy, and `pfpred verify`,
which checks every documented count and exits non-zero on drift. Phase 2 is in
progress: Task A clustering/splits have been run and validated on the supplied
Human Swiss-Prot FASTA (20,429 proteins, 13,070 clusters, zero cross-split
clusters); Task B (train-only vocabulary) remains.

---

## 8. Known limitations

- Track A uses contemporary labels and a clustered holdout; it is not a temporal
  or CAFA benchmark.
- Clustering at 30% identity and 80% bidirectional coverage still leaves
  sub-threshold homology; numbers remain
  optimistic relative to a truly novel protein.
- Track B follows CAFA6's temporal protocol; PK scores only terms not known at
  t0, and NK/LK are smaller.
- The frozen 2025-06 graph drops 0.06% of human annotations; irrelevant to Track
  A v1, relevant if Track A later uses IA F-max.
- ESM-2 is frozen; no task-specific fine-tuning in v1.
- Colab sessions are ephemeral; embeddings must be sharded and cached in Drive.
