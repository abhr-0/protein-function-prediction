# Data Sources and Provenance

All source data lives under `data/raw/` (gitignored) and is treated as
read-only.
There are no download scripts; the REST query and derivation rules below are
recorded for reproducibility.

---

## Track A — human Swiss-Prot

Source: UniProtKB REST API, current release **2026_04** (04-Nov-2026).

```
query:  reviewed:true AND organism_id:9606
url:    https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=txt
```

`format=txt` is used because it is the only format carrying both the sequence
and `DR GO` lines with aspect and UniProt evidence codes (`IDA`, `IPI`, `IEA`,
…). The TSV format exposes `go_id`, `go_f`, `go_p`, `go_c` but **no evidence**;
XML exposes ECO codes (`ECO:0000314`) instead of evidence codes.

The stream is 453 MB uncompressed. It is not kept. It is parsed once into:

| File | Size | Content |
|---|---|---|
| `data/raw/uniprot/human_proteins_2026_04.fasta.gz` | 6.4 MB | 20,429 sequences |
| `data/raw/uniprot/human_go_annotations_2026_04.tsv.gz` | 1.1 MB | 243,091 rows: `Entry`, `GO`, `Aspect`, `Evidence` |

Derivation rules (regeneration is a short one-off script, deliberately not
tracked):

1. Entry key = the first accession on the first `AC` line (the primary
   accession).
2. Sequence = the `SQ` block with all whitespace removed, headers dropped.
3. One row per `DR GO` line: `GO`, aspect letter, evidence code.
4. Evidence = the text after the last `; `, before `:` (`IDA:UniProtKB` →
   `IDA`).
5. No deduplication is needed: 243,091 rows are 243,091 unique protein–term
   pairs.
6. Output: UTF-8, tab-separated, gzip.

SHA-256 of the files on disk:

```
f66571d49bd900b2a8679dab4f5189d645b90bf729028ad12548b56d7b5c5d49  human_proteins_2026_04.fasta.gz
a7e993d9ab480952dbe7a6b4f87a580170577f9edc92e9edfe8ea30a6cae7cb0  human_go_annotations_2026_04.tsv.gz
```

Counts:

| Quantity | Value |
|---|---|
| Proteins | 20,429 |
| Residues | 11,417,710 (max 34,350) |
| Annotations (protein–term pairs, all evidence) | 243,091 |
| Terms | 17,317 |
| Proteins with at least one annotation | 19,375 |
| Annotations under `curated` policy | 216,272 |
| Proteins under `curated` policy | 18,630 |
| Terms under `curated` policy | 16,442 |
| Aspect counts (`curated`) | F 46,519 · P 96,961 · C 72,792 |

The `curated` policy excludes `IEA`, `ND`, and `NR`; it removes 26,819 IEA-only
pairs. The name avoids colliding with Track B's CAFA6 benchmark.

The UniProt web TSV `uniprotkb_human_AND_model_organism_9606_2026_11_18.tsv`
is **not** used: it has no evidence codes and no aspect column.

---

## Track B — CAFA6

Copied from the competition/evaluation bundle. The full test bundle is kept.

| File | Size | Purpose |
|---|---|---|
| `train/train_terms.tsv` | 10.7 MB | t0 labels |
| `train/train_sequences.fasta` | 53.3 MB | t0 sequences |
| `train/train_taxonomy.tsv` | 1.0 MB | species IDs |
| `test/testsuperset.fasta` | 101.1 MB | full test superset |
| `test/testsuperset-taxon-list.tsv` | 0.2 MB | species map |
| `eval/go-basic.obo` | 31.4 MB | frozen 2025-06 graph |
| `eval/IA.tsv` | 1.0 MB | official IA weights |
| `eval/test_sequences.fasta` | 102.9 MB | unchanged test sequences |
| `eval/May2026_groundtruth_{NK,LK,PK,PK_known}.tsv` | 5.6 MB | ground truth |
| `eval/May2026_groundtruth_terms_of_interest.txt` | 0.4 MB | terms of interest |
| `eval/README.md` | 5 KB | official protocol |
| `sample_submission.tsv` | 0.5 MB | submission format |

Counts:

| Quantity | Value |
|---|---|
| Train proteins / annotations / species | 82,404 / 537,027 / 1,381 |
| Ground-truth proteins | 8,023: NK 749, LK 485, PK 6,951 (162 in both LK and PK) |
| GT proteins present in t0 train | 7,219 of 8,023 |
| GT proteins present in Swiss-Prot 2026_04 | all 8,023; 2,232 are human |

`eval/README.md` states the same ground-truth total (8,023). `train_taxonomy.tsv`
has no header, so its 82,404 rows are all proteins.

**Leakage rule.** Track B labels must come from `train_terms.tsv` (t0). Do not
train Track B on the 2026_04 human files: 2,232 of the 8,023 GT proteins are
human and may carry post-t0 annotations. All 8,023 are present in the full
Swiss-Prot 2026_04 flat file.

---

## GO graph and IA

Track B uses CAFA6's frozen `eval/go-basic.obo` + `eval/IA.tsv`. No official IA
exists for 2026_04 labels. Track A v1 uses neither; it reports F1.

---

## Verification

`uv run pfpred verify` re-reads the raw files and compares 25 counts against the
constants in `pfpred.cli.EXPECTED`, exiting 1 on drift. Run it after any data
change. It takes about 5 seconds.

---

## Provenance rules

- `data/raw/` is read-only; never write derived data there.
- Reusable derived data (splits, vocabulary, embeddings) goes to `data/derived/`;
  run outputs go to `artifacts/runs/<run_id>/`.
- Record the REST release and query with any result.
- Hashes above identify the files on disk; regeneration is documented, not
  scripted.
