# Data Sources and Provenance

All large data files live under `data/` and are gitignored. Code, tests, and this
provenance note are tracked.

## Current snapshots

| File | Local path | Source | Size | Purpose |
|---|---|---|---|---|
| Swiss-Prot FASTA | `data/uniprot_sprot.fasta.gz` | [UniProt current release](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) | ~90 MB | Sequences |
| Swiss-Prot flat file | `data/uniprot_sprot.dat.gz` | [UniProt current release](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz) | ~700 MB | `DR GO` annotations with evidence codes |
| GO ontology | `data/go-basic.obo` | [GO current release](https://purl.obolibrary.org/obo/go/go-basic.obo) | ~32 MB | GO graph (`is_a`, `part_of`, `alt_id`) |
| CAFA6 benchmark | `data/cafa_6_protein_data/` | CAFA6 download | ~95 MB zip | Optional comparison mode, IA weights, submission format |

Detected release at the time of writing: **UniProt 2026_04 (04-Nov-2026)**,
576,073 Swiss-Prot entries. The local FASTA contains 576,073 records, so the
sequence snapshot matches that release.

## Downloading

```bash
# Default: flat file (GO annotations) + GO ontology
.venv/bin/python scripts/download_data.py

# Include the FASTA (re-download only if missing)
.venv/bin/python scripts/download_data.py --include dat fasta obo

# Force a refresh (may move to a newer release)
.venv/bin/python scripts/download_data.py --include dat obo --force

# Preview only
.venv/bin/python scripts/download_data.py --dry-run
```

The script writes a JSON manifest to `data/artifacts/data_manifest.json` with
each file's URL, size, SHA-256, detected release, and download timestamp.

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

UniProt GO annotations carry evidence codes. The planned training policy is:

- `cafa`: exclude `IEA` (automatic/electronic), `ND` (no data), and `NR`
  (not recorded). This approximates CAFA's training label policy.
- `full`: include IEA; still exclude `ND` and `NR`.

This is a working assumption from spot checks; `scripts/annotation_stats.py`
will report the evidence distribution so the policy can be confirmed before
training.

## Why the flat file instead of GOA GAF

`goa_uniprot_all.gaf.gz` is ~26 GB and must be filtered to reviewed accessions.
The Swiss-Prot flat file is ~700 MB and contains the same GO annotations with
evidence codes in `DR GO` lines. It also carries UniProt's own metadata and
secondary accessions, which helps join to the FASTA.

## CAFA6's role

CAFA6 files remain useful for:

- the IA weight table (`IA.tsv`),
- a controlled "CAFA labels vs UniProt labels" comparison,
- the submission format for `testsuperset.fasta`.

They are not the default training label source. CAFA6 test labels are withheld,
and `testsuperset.fasta` contains all training proteins, so local CAFA6 test
scoring is not possible from the released files.
