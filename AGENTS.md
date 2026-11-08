# AGENT.md — Operating Guidelines for AI Agents

> **Project:** Protein Function Predictor (UniProtKB/Swiss-Prot & CAFA Benchmark)  
> **Primary References:** [`protein_function_predictor_knowledge_and_plan.md`](./protein_function_predictor_knowledge_and_plan.md), UniProtKB Swiss-Prot  
> **Status:** Phase 1 (Data Prep) & Phase 2 (Classical ML Baseline)  

---

## 1. Project Overview & Objective

The objective of this project is to build an AI-based system that predicts biological functions (Gene Ontology terms: Molecular Function `MF`, Biological Process `BP`, Cellular Component `CC`) directly from protein amino-acid sequences.

### Primary Dataset: UniProtKB / Swiss-Prot
The project operates directly on the **UniProtKB/Swiss-Prot** dataset ([`uniprot_sprot.fasta.gz`](./uniprot_sprot.fasta.gz)):
- **Scope:** 575,503 high-quality, manually reviewed, non-redundant protein sequences.
- **FASTA Header Schema:** `>db|Accession|EntryName Description OS=Organism OX=TaxonID [GN=Gene] PE=Evidence SV=Version`
  - Example: `>sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 OX=654924 GN=FV3-001R PE=4 SV=1`
  - Primary Key (`EntryID`): The UniProt Accession ID (e.g., `Q6GZX4`).
- **Annotations:** Gene Ontology (GO) annotations mapped via UniProt accession IDs (using `cafa_6_protein_data/Train/train_terms.tsv` which covers 82,396 Swiss-Prot proteins with 537k+ validated GO terms, with option to extend to full UniProt GOA).

### Modeling Strategy
1. **Stage 1 (Current Focus):** Classical Machine Learning baseline using hand-engineered sequence features (amino acid composition, length, physicochemical properties, k-mers) predicting multi-label GO terms. Focus on **Molecular Function (MF)** first.
2. **Stage 2 (Future Focus):** Pretrained protein language model embeddings (e.g., ESM-2) extracted from UniProt sequences.
3. **Core Research Question:** Comparative evaluation of classical hand-crafted features vs. learned PLM embeddings on Swiss-Prot sequences using CAFA-style metrics (Information Accretion weighted F-max).

---

## 2. Hardware & Environment Constraints

- **Dataset Scale Management:**
  - `uniprot_sprot.fasta.gz` contains 575k sequences (~90 MB compressed, ~300 MB uncompressed).
  - **Streaming from Gzip:** Never decompress the full file to disk; stream directly using `gzip.open('uniprot_sprot.fasta.gz', 'rt')` and `Bio.SeqIO.parse`.
  - **Dev Mode / Subsampling:** For rapid prototyping and testing, always support a `dev_mode` or `max_records` flag (e.g., ~5,000 proteins, top 100–300 most frequent MF GO terms).
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

# Run a python script in the environment
.venv/bin/python src/data/loader.py
```

---

## 4. Repository Layout & Conventions

```text
protein-function-prediction/
├── data/                      # Unified data directory (ignored by git except .gitkeep)
│   ├── uniprot_sprot.fasta.gz # Primary UniProtKB/Swiss-Prot sequence database (compressed)
│   ├── cafa_6_protein_data/   # CAFA benchmark annotations & evaluation files
│   │   ├── Train/
│   │   │   ├── train_terms.tsv    # UniProt EntryID -> GO term annotations
│   │   │   ├── train_taxonomy.tsv # Taxon mappings
│   │   │   └── go-basic.obo       # Gene Ontology graph hierarchy
│   │   ├── IA.tsv                 # Information Accretion weights for CAFA eval
│   │   └── sample_submission.tsv
│   └── artifacts/             # Saved models, precomputed feature matrices, metrics
├── notebooks/                 # Exploratory data analysis & visualizations
│   └── 1_preprocessing.ipynb
├── src/                       # Modular production pipeline
│   ├── data/                  # UniProt FASTA streaming, header parser, annotation joiner
│   ├── features/              # Feature extractors (AAC, physicochemical, k-mers)
│   ├── models/                # Baseline classifiers (Logistic Regression, LightGBM, etc.)
│   ├── evaluation/            # Multi-label metrics, CAFA F-max, IA weighting
│   └── utils/                 # Path handling, logging, caching
├── tests/                     # Automated unit and integration tests (Pytest)
├── pyproject.toml             # Project dependencies & metadata
├── uv.lock                    # Dependency lockfile
├── flake.nix                  # Nix dev environment
├── protein_function_predictor_knowledge_and_plan.md # Reference doc
└── AGENT.md                   # This operational manual
```

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
5. **Memory Discipline:**
   - Never load 575k dense feature vectors into memory simultaneously. Use chunked processing or sparse matrices.
6. **Git Hygiene:**
   - Keep `uniprot_sprot.fasta.gz` and cache files out of git commits. Track code and configurations only.
