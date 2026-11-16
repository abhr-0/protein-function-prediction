# AGENTS.md — Operating Guidelines for AI Agents

> **Project:** Protein Function Predictor (UniProtKB/Swiss-Prot & CAFA Benchmark)  
> **Primary References:** UniProtKB/Swiss-Prot release notes; `ROADMAP.md` for phase status and measured facts  
> **Status:** Phase 1 (UniProt labels + leakage-safe splits) complete; Phase 2 (GO propagation + IA-weighted F-max) pending  

---

## 1. Project Overview & Objective

The objective of this project is to build an AI-based system that predicts biological functions (Gene Ontology terms: Molecular Function `MF`, Biological Process `BP`, Cellular Component `CC`) directly from protein amino-acid sequences.

### Primary Dataset: UniProtKB / Swiss-Prot
The project operates directly on the **UniProtKB/Swiss-Prot** dataset (`data/uniprot/uniprot_sprot.fasta.gz`):
- **Scope:** 576,073 high-quality, manually reviewed, non-redundant protein sequences (UniProt release 2026_04).
- **FASTA Header Schema:** `>db|Accession|EntryName Description OS=Organism OX=TaxonID [GN=Gene] PE=Evidence SV=Version`
  - Example: `>sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 OX=654924 GN=FV3-001R PE=4 SV=1`
  - Primary Key (`EntryID`): The UniProt Accession ID (e.g., `Q6GZX4`).
- **Annotations (primary):** Gene Ontology terms parsed from the UniProt flat file `data/uniprot/uniprot_sprot.dat.gz` `DR GO` lines, from the **same release** as the FASTA (2026_04, 3,100,768 raw GO lines). Evidence policy `cafa` excludes `IEA`/`ND`/`NR`; `full` keeps `IEA`.
- **Benchmark labels (secondary):** `cafa_6_protein_data/Train/train_terms.tsv` covers only 82,395 of 576,073 Swiss-Prot proteins (**14.3%**) and its test labels are withheld, so it is kept only as an optional comparison/submission mode — never as the default training source.

### Modeling Strategy
1. **Stage 1 (Current Focus):** Classical Machine Learning baseline using hand-engineered sequence features (amino acid composition, length, physicochemical properties, k-mers) predicting multi-label GO terms. Focus on **Molecular Function (MF)** first.
2. **Stage 2 (Future Focus):** Pretrained protein language model embeddings (e.g., ESM-2) extracted from UniProt sequences.
3. **Core Research Question:** Comparative evaluation of classical hand-crafted features vs. learned PLM embeddings on Swiss-Prot sequences using CAFA-style metrics (Information Accretion weighted F-max).

---

## 2. Hardware & Environment Constraints

- **Dataset Scale Management:**
  - `uniprot_sprot.fasta.gz` contains 576k sequences (~90 MB compressed).
  - `uniprot_sprot.dat.gz` (labels) is ~668 MB compressed and must be streamed, never decompressed to disk.
  - **Streaming from Gzip:** Never decompress the full file to disk; stream directly using `gzip.open('uniprot_sprot.fasta.gz', 'rt')` and `Bio.SeqIO.parse`.
  - **Dev Mode / Subsampling:** For rapid prototyping use a size limit (default `--n-proteins 5000`, MF). Define the label space by training frequency (`--min-term-count`, default 20) rather than a fixed top-K cap: only ~52 MF terms reach 20 training proteins at 5k, so a fixed cap of 100–300 either cannot be reached or admits singleton terms.
  - **Sparse Representation:** Always store multi-label target matrices and high-dimensional feature matrices as `scipy.sparse.csr_matrix` to keep memory consumption under 2 GB.

---

## 3. Package & Environment Management

The repository uses **Nix** and **`uv`** for dependency and environment management.

- **Python Virtualenv:** `.venv` (Python 3.13+)
- **Executable path:** `.venv/bin/python`, `.venv/bin/pytest`
- **Package configuration:** `pyproject.toml` and `uv.lock`
- **Installed Dependencies:** `scikit-learn`, `scipy`, `biopython`, `pandas`, `pytest`, `jupyter`
- **Nix Flake:** `flake.nix` provides the reproducible shell with dynamic library paths (`LD_LIBRARY_PATH`).

### Key Commands
```bash
# Run unit tests
.venv/bin/pytest

# Add dependencies with uv (inside nix develop shell)
uv add <dependency>

# Run a pipeline script in the environment
.venv/bin/python scripts/run_baseline_subset.py --label-source cafa --n-proteins 200

# Report label coverage for an evidence policy
.venv/bin/python scripts/annotation_stats.py --aspect F --evidence-policy cafa
```

---

## 4. Repository Layout & Conventions

```text
protein-function-prediction/
├── data/                      # Unified data directory (gitignored)
│   ├── uniprot/               # UniProt-published data (one source, one release)
│   │   ├── uniprot_sprot.fasta.gz # Swiss-Prot sequences (2026_04)
│   │   └── uniprot_sprot.dat.gz   # Swiss-Prot flat file: DR GO labels + evidence
│   ├── go/                    # GO Consortium data (own release cadence)
│   │   └── go-basic.obo       # GO graph for ancestor propagation (Phase 2)
│   ├── cafa_6_protein_data/   # Optional comparison mode only
│   │   ├── Train/
│   │   │   ├── train_terms.tsv    # CAFA6 labels (14.3% Swiss-Prot coverage)
│   │   │   ├── train_sequences.fasta
│   │   │   ├── train_taxonomy.tsv
│   │   │   └── go-basic.obo       # CAFA6's own OBO snapshot (do not mix with the fresh one)
│   │   ├── IA.tsv                 # Information Accretion weights for CAFA eval
│   │   └── sample_submission.tsv
│   └── artifacts/             # Models, metrics JSON, data_manifest.json, annotation cache
├── docs/
│   └── data_sources.md        # URLs, release, evidence policy, provenance rules
├── notebooks/                 # Exploratory data analysis & visualizations
├── scripts/                   # Thin CLI wrappers over src/
│   ├── download_data.py       # Fetch external data + write provenance manifest
│   ├── annotation_stats.py    # Coverage / evidence distribution report
│   └── run_baseline_subset.py # End-to-end baseline run
├── src/                       # Modular production pipeline
│   ├── data/
│   │   ├── uniprot.py         # FASTA streaming, header parsing, alignment
│   │   ├── annotations.py     # UniProt DR GO parsing, evidence policy, cache
│   │   ├── targets.py         # Train-only vocabulary, apply_vocabulary, coverage
│   │   └── prepare.py         # split -> vocabulary -> filter (the leakage guard)
│   ├── features/              # Feature extractors (AAC, physicochemical, k-mers)
│   ├── models/                # Baseline classifiers (Logistic Regression, LightGBM, etc.)
│   ├── evaluation/            # Multi-label metrics, CAFA F-max, IA weighting
│   └── utils/                 # Path handling, download, provenance
├── tests/                     # Automated unit and integration tests (Pytest, 104 tests)
├── pyproject.toml             # Project dependencies & metadata
├── uv.lock                    # Dependency lockfile
├── flake.nix                  # Nix dev environment
├── ROADMAP.md                 # Phase status, locked decisions, measured facts
└── AGENTS.md                  # This operational manual
```

### Planning documents

- **`ROADMAP.md` (tracked)** answers *what is the project's direction and which
  milestones are closed*. It changes only at milestone boundaries and is
  authoritative for phase status.
- **`PLAN.md` (untracked, optional)** is a local working plan for the current
  phase: commit-level checklists, open questions, scratch measurements. **If
  `PLAN.md` exists in the project root, read it before starting work** — it
  usually describes the immediate task in more detail than `ROADMAP.md`.
  Because it is untracked it may be stale, absent on a fresh clone, or
  partially finished: it is subordinate to `ROADMAP.md` and this file. If they
  disagree, follow `ROADMAP.md` and say so rather than silently picking one.

---

## 5. Development Principles & Agent Rules

When writing code or planning changes in this repo, agents must adhere to the following rules:

1. **Native UniProt Stream Ingestion:**
   - Always parse UniProt records directly from `uniprot_sprot.fasta.gz` using generators/streaming.
   - Parse headers to extract: `EntryID` (Accession), `EntryName` (e.g. `001R_FRG3G`), `Organism` (`OS=`), `TaxonID` (`OX=`), and `Sequence`.
2. **Modular Code Over Monolithic Notebooks:**
   - All data parsing, feature extraction, training, and metrics must live in `src/`.
   - Notebooks in `notebooks/` are strictly for visualization, EDA, and inspection.
3. **Verify with Tests (`pytest`):**
   - Every module in `src/` must have a corresponding test in `tests/` tested against small synthetic data or a tiny subset (e.g., 5 sample UniProt records).
   - Run `.venv/bin/pytest` and verify pass before reporting task completion.
4. **Prevent Data Leakage:**
   - Train/Val/Test splits must be made at the protein `EntryID` level before any feature engineering, vocabulary building, or scaling.
   - Term vocabulary is chosen on the **training split only** (`select_target_terms`), then applied to every split (`apply_vocabulary`). Never filter terms inside `load_annotations()` or `build_annotated_dataset()`; both are intentionally vocabulary-free.
   - The zero-positive filter (`apply_vocabulary`) must run after the split and after vocabulary selection, and its drop counts must be reported, not hidden.
   - `prepare_splits()` owns that ordering end to end; call it instead of re-implementing split/vocabulary/filter logic in scripts.
5. **Memory Discipline:**
   - Never load 576k dense feature vectors into memory simultaneously. Use chunked processing or sparse matrices.
6. **Git Hygiene:**
   - Keep `data/uniprot/`, `data/go/`, `data/cafa_6_protein_data/`, `data/artifacts/`, and cache files out of git commits. Track code, docs, and configurations only.
