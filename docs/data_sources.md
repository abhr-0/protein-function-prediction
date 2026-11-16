# Data Sources and Provenance

All large data files live under `data/` and are gitignored. Code, tests, and this
provenance note are tracked.

## Current snapshots

| File | Local path | Source | Size | Purpose |
|---|---|---|---|---|
| Swiss-Prot FASTA | `data/uniprot/uniprot_sprot.fasta.gz` | [UniProt current release](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) | ~90 MB | Sequences |
| Swiss-Prot flat file | `data/uniprot/uniprot_sprot.dat.gz` | [UniProt current release](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz) | 668 MB | `DR GO` annotations with evidence codes |
| GO ontology | `data/go/go-basic.obo` | [GO current release](https://purl.obolibrary.org/obo/go/go-basic.obo) | 31 MB | GO graph (`is_a`, `part_of`, `alt_id`) — see the note on two OBO copies below |
| CAFA6 train/test | `data/cafa6/{train,test}/` | [Kaggle competition](https://www.kaggle.com/competitions/cafa-6-protein-function-prediction) (login + rules) | ~161 MB | Optional comparison mode, submission targets |
| CAFA6 evaluation | `data/cafa6/eval/` | [figshare, DOI 10.6084/m9.figshare.32861969](https://doi.org/10.6084/m9.figshare.32861969) (CC BY 4.0) | ~135 MB | Official `IA.tsv`, evaluation graph, May 2026 ground truth |

Release: **UniProt 2026_04 (04-Nov-2026)**, 576,073 Swiss-Prot entries.

Verified on 2026-11-16, so sequences and labels are provably one snapshot:

| Check | Result |
|---|---|
| FASTA records | 576,073 |
| Flat-file entries (`^ID`) | 576,073 |
| Raw `DR GO` lines | 3,100,768 |
| Flat-file md5 vs upstream `RELEASE.metalink` | matches (`05779136a543d28d8d0ae9488d313a5e`) |

The manifest stores SHA-256 for each file, but UniProt's metalink publishes
**md5 only** for the flat file, so the manual check above is md5-strength. The
script does not yet verify against upstream hashes (listed under gaps).

## Downloading

```bash
# Default: labels + sequences + GO ontology
.venv/bin/python scripts/download_data.py

# Add the official CAFA6 evaluation bundle (figshare, md5-verified)
.venv/bin/python scripts/download_data.py --include cafa6

# Everything the script can fetch
.venv/bin/python scripts/download_data.py --include all

# Force a refresh (may move to a newer release)
.venv/bin/python scripts/download_data.py --force
```

The script writes `data/artifacts/data_manifest.json` with each file's URL,
size, SHA-256, md5 (when upstream publishes one), publisher label, release, and
download timestamp. Manifest keys are publisher-prefixed
(`uniprot/…`, `go/…`, `cafa6-figshare/…`) because two different files are both
named `go-basic.obo`.

Releases are recorded **per publisher**: UniProt files carry the UniProt
release, the ontology carries its own header version, and the CAFA6 evaluation
bundle carries `CAFA6 May2026`. Files whose upstream md5 is known are verified
after transfer; a mismatch deletes the partial file and fails the run.

## Provenance rule

**Sequences and labels must come from the same snapshot.** If the FASTA and the
flat file are from different releases, re-download the cheaper FASTA so they
match:

```bash
.venv/bin/python scripts/download_data.py --include fasta --force
```

Do not mix CAFA6 labels (an older snapshot) with current Swiss-Prot sequences
unless the mismatch is deliberate and documented. For the CAFA6 comparison
mode, use CAFA6 `train_sequences.fasta` to guarantee exact snapshot consistency.

## Evidence policy

UniProt GO annotations carry evidence codes. The training policy is:

- `cafa`: exclude `IEA` (automatic/electronic), `ND` (no data), and `NR`
  (not recorded). This approximates CAFA's training label policy.
- `full`: include IEA; still exclude `ND` and `NR`.

Under the `cafa` policy the file yields **1,292,302 annotations over 183,145
proteins and 27,559 GO terms** (`.venv/bin/python scripts/annotation_stats.py`,
2026-11-16). Per aspect:

| Aspect | Annotations | Proteins | Terms |
|---|---|---|---|
| F (MF) | 293,489 | 135,237 | 7,649 |
| P (BP) | 577,377 | 135,800 | 17,004 |
| C (CC) | 421,436 | 153,700 | 2,906 |

**Still open:** the claim that CAFA6's own `train_terms.tsv` excludes `IEA` was
inferred from spot checks (only ~2.5% ancestor slots, no `IEA`-only terms), not
verified against CAFA's published policy. The coverages above confirm the
policy leaves plenty of labels; they do not confirm it matches CAFA.

## Two copies of `go-basic.obo`

`data/go/go-basic.obo` (fresh, used for propagation) and
`data/cafa6/eval/go-basic.obo` (CAFA6's frozen evaluation graph, GO release
2025-06-01) differ in
vintage, and therefore in `alt_id` and obsolete-term handling. Rule: propagate
with the **fresh** copy; use CAFA6's copy only when reproducing CAFA6's own
setup, and never mix the two in one run.

## Why the flat file instead of GOA GAF

`goa_uniprot_all.gaf.gz` is ~26 GB and must be filtered to reviewed accessions.
The Swiss-Prot flat file is ~700 MB and contains the same GO annotations with
evidence codes in `DR GO` lines. It also carries UniProt's own metadata and
secondary accessions, which helps join to the FASTA.

## CAFA6's role

CAFA6 ran on Kaggle from 2025-10-15 to 2026-06-01. The training bundle covers
82,395 of 576,073 Swiss-Prot proteins (14.3%), so it is never the default
training source. Two separate things are provided:

**Competition bundle (Kaggle, login required)** — `train_terms.tsv`,
`train_sequences.fasta`, `train_taxonomy.tsv`, `testsuperset.fasta`,
`testsuperset-taxon-list.tsv`, `sample_submission.tsv`. Verified authentic
against the competition's published breakdown (537,027 annotations:
BP 250,805 / CC 157,770 / MF 128,452).

**Evaluation bundle (figshare, open)** — everything needed to score locally:

| File | Purpose |
|---|---|
| `IA.tsv` | information accretion weights, matching the evaluation graph |
| `go-basic.obo` | frozen evaluation graph, GO release **2025-06-01** |
| `test_sequences.fasta` | 224,121 proteins whose sequences were unchanged between t0 and t1 |
| `May2026_groundtruth_{NK,LK,PK,PK_known}.tsv` | ground truth per evaluation regime |
| `May2026_groundtruth_terms_of_interest.txt` | terms valid in both the frozen and t1 graph |

**Ground truth is public.** During the competition these labels were secret;
they were released on 2026-07-01, so CAFA6 can now be scored locally instead of
only submitted. Regime sizes: NK 749, LK 485, PK 6,951 proteins (8,023 unique).

**Official protocol** (from the bundle's README):

```bash
cafaeval go-basic.obo predictions/ May2026_groundtruth_PK.tsv \
  -ia IA.tsv -known May2026_groundtruth_PK_known.tsv \
  -no_orphans -norm cafa -prop fill \
  -toi May2026_groundtruth_terms_of_interest.txt \
  -max_terms 500 -th_step 0.001
```

Only PK uses `-known`. The ground truth deliberately excludes
UniProtKB:SubCell and RHEA annotations.

## Known gaps in this provenance

- The CAFA6 **training/test** bundle still requires a Kaggle login and rules
  acceptance, so `scripts/download_data.py` cannot fetch it. Keep the local copy
  and check it against the counts above.
- Downloads record upstream md5 for sources that publish one (CAFA6 figshare) but
  UniProt's `RELEASE.metalink` md5 is not yet wired into the script, so the
  UniProt files remain manually verified.
- The CAFA6 evaluation graph (2025-06-01) and the fresh GO graph are different
  vintages; use each only for its own purpose.
