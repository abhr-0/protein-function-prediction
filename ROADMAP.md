# ROADMAP.md — Project Direction and Milestone Status

> Durable view of where this project is going and which milestones are closed.
> Last updated: 2026-11-16.
>
> **Division of labour:** this file holds milestone-level status, locked decisions,
> and measured facts — it should change only at milestone boundaries.
> Commit-level checklists, open questions, and day-to-day detail live in the
> working plan (`PLAN.md`); operational rules for agents live in `AGENTS.md`.

---

## Status at a glance

| Phase | Scope | Status |
|---|---|---|
| 0. Foundation | Environment, ingestion, features, baseline model, basic metrics | **Complete** |
| 1. Label correctness | Swiss-Prot labels from UniProt `DR GO`, leakage-safe splits | **Complete (2026-11-16)** |
| 2. CAFA-style evaluation | GO propagation + IA-weighted protein-centric F-max | **Not started** |
| 3. Extra signals | CC aspect, localization ablation, k-NN, PPI, CAFA submission | **Not started** |

Test suite: **104 tests passing** (`.venv/bin/pytest`).

---

## Phase 0 — Foundation (complete)

- [x] Nix + `uv` environment; `scikit-learn`, `scipy`, `biopython`, `pandas`, `pytest`.
- [x] Package skeleton: `src/data`, `src/features`, `src/models`, `src/evaluation`, `src/utils`.
- [x] Streaming Swiss-Prot ingestion (`stream_uniprot_fasta`, `parse_uniprot_header`).
- [x] Features: 20-dim AAC, 6 physicochemical descriptors, sparse 2-mers, unified pipeline (~426 columns).
- [x] Baseline classifier: sparse target matrix, `EntryID`-level splits, Binary Relevance logistic regression.
- [x] Basic metrics: micro/macro P-R-F1, deterministic validation threshold selection.

## Phase 1 — Label correctness (complete)

The defect being fixed: training on CAFA6 `train_terms.tsv` covered only 82,395 of
576,073 Swiss-Prot proteins (14.3%), and its vocabulary was chosen before the split.

- [x] Download helper with SHA-256 + provenance manifest; UniProt 2026_04 FASTA and
      flat file verified as the **same snapshot** (576,073 entries each; local `.dat.gz`
      md5 matches the upstream `RELEASE.metalink` value).
- [x] Streaming `DR GO` parser with evidence policies (`cafa` excludes `IEA`/`ND`/`NR`;
      `full` keeps `IEA`) and a gzip annotation cache.
- [x] Annotation coverage/evidence report (`scripts/annotation_stats.py`).
- [x] Source-agnostic alignment builder (`build_annotated_dataset`) with seeded uniform
      protein sampling.
- [x] Train-only vocabulary selection, single zero-positive filter
      (`select_target_terms` → `apply_vocabulary`), coverage reporting.
- [x] Runner wiring: `--label-source {uniprot,cafa}`, 3-way split, threshold tuned on
      validation and reported on test, provenance JSON.
- [x] **Record the Phase 1 sanity run** (2026-11-16), UniProt labels, 5k proteins, MF
      (the current defaults reproduce it exactly):
      `.venv/bin/python scripts/run_baseline_subset.py --label-source uniprot \
      --n-proteins 5000 --aspect F`
      → 3,500 / 500 / 1,000 train/val/test, 52-term vocabulary, test annotation coverage
      33.5%, validation-tuned threshold 0.11, **test micro F1 0.209** (P 0.189 / R 0.235),
      macro F1 0.095. Plumbing check only — **not** a CAFA result.

Acceptance criteria for closing Phase 1:

- [x] `.venv/bin/pytest` passes.
- [x] Old `--label-source cafa` path still runs end-to-end.
- [x] `--label-source uniprot --n-proteins 5000 --aspect F` runs end-to-end.
- [x] Coverage report shows the expected label counts: 1,292,302 curated annotations over
      183,145 proteins and 27,559 terms (`cafa` policy); MF alone is 293,489 annotations
      over 135,237 proteins.
- [x] Every selected term has at least one positive in the training matrix (guaranteed by
      construction: the vocabulary is derived from training counts with `min_count >= 1`).
- [x] Comparison mode still works on current code: `--label-source cafa --n-proteins 2000
      --min-term-count 5 --top-terms 100` completes with a 76-term vocabulary and test micro
      F1 0.629 (a different label space, so not comparable to the UniProt number).
- [x] No test protein contributes to the vocabulary or to feature fitting — enforced by
      `prepare_splits()` (split -> train-only vocabulary -> restrict every split) and
      asserted by `tests/test_prepare.py::test_vocabulary_comes_from_training_terms_only`;
      the feature half is covered by
      `test_kmer_vocabulary_is_fitted_on_training_sequences_only`. A runner-level
      integration test (`tests/test_runner_integration.py`) exercises the wiring on
      synthetic data.

### Bugs found by the sanity run (both fixed)

- GO annotations were keyed on the **last** `AC` line of an entry instead of the primary
  accession. Merged entries wrap `AC` across several lines, so annotations were attributed
  to a secondary accession absent from the FASTA: 63 of the first 5,000 sampled proteins
  were silently dropped, and annotation rows collided during de-duplication.
- `build_annotated_dataset` crashed on the UniProt cache's categorical dtypes
  (`TypeError: unhashable type: 'list'`), because a list-valued aggregation cannot be cast
  back to a categorical column.

## Phase 2 — CAFA-style evaluation (not started)

- [ ] GO OBO parser + ancestor propagation (`is_a`, `part_of`, `alt_id`, obsolete terms);
      stdlib only, no new dependency. Evaluation-side first; training-side propagation is
      a later ablation.
- [ ] IA-weighted protein-centric F-max: propagate ground truth as a set, propagate
      predictions by max score over descendants, sweep thresholds, report the maximum.
- [ ] Wire into the runner with per-aspect (F/P/C) reporting, label-space coverage, and
      provenance in the results JSON. IA weights come from CAFA `IA.tsv`.

## Phase 3 — Extra signals (after Phases 1–2 are stable)

- [ ] Cellular Component as a third prediction aspect (73% of CAFA-annotated proteins
      carry CC terms, so CC is not sparse). Per-aspect protein sets; do not intersect
      proteins across aspects.
- [ ] Predicted localization as an input feature — optional ablation only; never use
      curated CC labels as inputs.
- [ ] Sequence-similarity k-NN baseline (DIAMOND/BLAST test-vs-train, aggregate
      neighbour GO terms). Sequence-only, so it stays comparable to the feature model.
- [ ] CAFA6 comparison mode / submission writer (`sample_submission.tsv` format). No
      local scoring: CAFA6 test labels are withheld.
- [ ] PPI network-assisted track — **deferred**. Separate research question requiring
      external data, ID mapping, and strict leakage controls; must not be mixed with the
      sequence-only comparison.

---

## Measured facts

| Quantity | Value |
|---|---|
| Swiss-Prot entries (release 2026_04) | 576,073 |
| Raw `DR GO` lines in `uniprot_sprot.dat.gz` | 3,100,768 |
| CAFA6 `train_terms.tsv` proteins | 82,404 (82,395 present in Swiss-Prot) |
| Swiss-Prot covered by CAFA6 labels | **14.3%** |
| Proteins dropped by the old CAFA-only pipeline | ~493,678 |
| CAFA6 test labels | withheld (no local CAFA6 test scoring possible) |

## Locked decisions

| Decision | Choice |
|---|---|
| Training labels | UniProt `DR GO` from `uniprot_sprot.dat.gz`, same release as the FASTA |
| Evidence policy | `cafa` = exclude `IEA`/`ND`/`NR`; `full` = include `IEA`, still drop `ND`/`NR` |
| Benchmark | Local `EntryID`-level Swiss-Prot holdout with CAFA IA-weighted F-max |
| CAFA6 role | Optional comparison / submission mode only |
| Vocabulary | Selected on the training split only, then applied to every split |
| Propagation | Evaluation-side by default; training-side is an ablation |
| Model/features during the label fix | Unchanged, so metric changes are attributable to the data |
| GO graph for propagation | `data/go-basic.obo` (fresh). CAFA6's older copy is only for reproducing CAFA6's setup; never mix the two in one run |

## Known limitations

- The random holdout is optimistic: homologs can appear on both sides of the split.
  Sequence-cluster or temporal splits are future work.
- Contemporary labels plus a random holdout is not CAFA's future-annotation setting, so
  absolute numbers will not be comparable to published CAFA leaderboards.
- `scripts/download_data.py` does not yet fetch `cafa_6_protein_data/`, so CAFA mode is
  not reproducible from the script alone.
- **MF labels are extremely top-heavy**: only 52 MF terms occur in >=20 of 3,500 training
  proteins. A small vocabulary caps coverage (52 terms -> 33.5% of test proteins), while a
  wide one (287 terms at `min_count=5`) drops sparsity to 99.4% and makes the classifier's
  probabilities too low to threshold.
- **Global micro-F1 threshold selection is fragile** with many rare labels: it previously
  chose 0.00 ("predict everything", precision 0.006) on a 290-term vocabulary. The grid now
  excludes 0.0 and 1.0, and a threshold at the search floor is flagged in the console and in
  the metrics JSON (`at_search_floor`). This does not make the model better — the same run
  now reports F1 0.001 instead of an inflated 0.012 — it stops the metric rewarding an
  uninformative all-positive rule. Per-class thresholds and the Phase 2 IA-weighted F-max
  remain the real remedies.
- Artifact names now encode the vocabulary policy
  (`baseline_metrics_<source>_<n>_<aspect>_mc<min_count>_k<cap>.json`), so runs with
  different settings no longer overwrite each other. The seed is still not in the name.
