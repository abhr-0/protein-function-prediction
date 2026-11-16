# Data Sources and Provenance

All large data files live under `data/` and are gitignored. Code, tests, and this
provenance note are tracked.

## Current snapshots

| File | Local path | Source | Size | Purpose |
|---|---|---|---|---|
| Swiss-Prot FASTA | `data/uniprot/uniprot_sprot.fasta.gz` | [UniProt current release](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) | ~90 MB | Sequences |
| Swiss-Prot flat file | `data/uniprot/uniprot_sprot.dat.gz` | [UniProt current release](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz) | 668 MB | `DR GO` annotations with evidence codes |
| GO ontology | `data/go/go-basic.obo` | [GO current release](https://purl.obolibrary.org/obo/go/go-basic.obo) | 31 MB | GO graph (`is_a`, `part_of`, `alt_id`) — see the note on two OBO copies below |
| CAFA6 benchmark | `data/cafa_6_protein_data/` | CAFA6 download (**URL not yet recorded — see gaps**) | ~95 MB zip | Optional comparison mode, IA weights, submission format |

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
`data/cafa_6_protein_data/Train/go-basic.obo` (CAFA6's older snapshot) differ in
vintage, and therefore in `alt_id` and obsolete-term handling. Rule: propagate
with the **fresh** copy; use CAFA6's copy only when reproducing CAFA6's own
setup, and never mix the two in one run.

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
scoring is not possible from the released files. CAFA6 can therefore only ever be
a submission target, never a locally scored benchmark.

## Known gaps in this provenance

- `scripts/download_data.py` does not fetch `cafa_6_protein_data/` at all, so the
  IA weights that Phase 2's metric depends on are not reproducible from the
  script. The CAFA6 download URL is not recorded anywhere yet.
- `fasta` is not in the default `--include` list, so a fresh clone gets labels and
  ontology but no sequences.
- The manifest records the **UniProt** release on the `go-basic.obo` entry, but GO
  has its own release cadence. That field is misleading for the one file whose
  version decides propagation behaviour.
- Downloads are hashed after transfer but never checked against upstream hashes
  (`RELEASE.metalink` for UniProt).
