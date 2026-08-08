# Output formats

This document describes CHIMERA bundle schema version
`urn:chimera:benchmark-bundle:1` and split schema
`urn:chimera:split-manifest:1`. Treat schema URNs and column names—not prose
layout or file ordering—as the machine interface. Validate before consuming:

```console
chimera validate PATH/TO/BUNDLE
chimera validate PATH/TO/BUNDLE --json
```

## Bundle tree

A full `chimera suite` writes:

```text
bundle/
├── .chimera-bundle
├── manifest.json
├── resolved-config.json
├── execution.json
├── references.tsv
├── sequences.tsv
├── source-sequences.fasta.gz
├── schemas/
├── excluded.tsv
├── REPORT.md
├── checksums.sha256
├── 2a_random_fragment/
│   ├── assignments.tsv
│   ├── excluded.tsv
│   ├── split.json
│   ├── train.fasta.gz
│   ├── train.truth.tsv.gz
│   ├── test.fasta.gz
│   └── test.truth.tsv.gz
├── 2b_genome_holdout/
│   └── (the same seven protocol files)
├── 2c_similarity_filtered/
│   ├── assignments.tsv
│   ├── excluded.tsv
│   ├── split.json
│   ├── train.fasta.gz
│   ├── train.truth.tsv.gz
│   ├── test.fasta.gz
│   ├── test.truth.tsv.gz
│   ├── candidate_test.fasta.gz
│   ├── candidate_test.truth.tsv.gz
│   ├── external-similarity.tsv  (only when an external table was supplied)
│   └── test_strata/
│       ├── high_similarity.fasta.gz
│       ├── high_similarity.truth.tsv.gz
│       ├── moderate_similarity.fasta.gz
│       ├── moderate_similarity.truth.tsv.gz
│       ├── low_similarity.fasta.gz
│       ├── low_similarity.truth.tsv.gz
│       ├── distant_detectable.fasta.gz
│       ├── distant_detectable.truth.tsv.gz
│       ├── no_detectable_match.fasta.gz
│       └── no_detectable_match.truth.tsv.gz
├── 2d_temporal_holdout/
│   └── (the same seven protocol files)
└── 2e_taxonomic_holdout/
    └── (the same seven protocol files)
```

If `generate` selects fewer protocols, only their stable directories are
present. Empty exclusion tables and empty similarity strata are still written
with headers. `external-similarity.tsv` is present only when the resolved
configuration supplies an external similarity table; it is the verbatim
evidence snapshot used for Test 2C.

Version 1 has a closed filesystem layout. `chimera validate` rejects unknown
files or directories anywhere in the bundle, including entries added to
`checksums.sha256`. The conditional external-similarity snapshot above is the
only protocol-specific optional file.

## Serialization conventions

- Text and tables are UTF-8 with LF line endings.
- TSV uses a header row and RFC-4180-style quoting with tab as delimiter.
- Dates use ISO `YYYY-MM-DD`; timestamps use ISO 8601 with an explicit UTC
  offset.
- Missing optional TSV values are empty fields, never textual `NA`, `null`, or
  `None`.
- Similarity and coverage values are fractions on `[0,1]`, serialized with the
  shortest decimal representation that round-trips to the same IEEE-754 value.
  Scientific evidence is never quantized for presentation.
- JSON is UTF-8, indented, key-sorted, and rejects NaN/infinity.
- FASTA and compressed truth are gzip streams with `mtime=0` and no embedded
  filename. FASTA sequence lines are wrapped at 80 characters.
- Generated FASTA headers contain only the opaque `sequence_id`. Join to truth
  instead of parsing an identifier.
- Fragment coordinates are 0-based and half-open on the original forward
  source contig. Circular intervals are unwrapped and may have an end greater
  than source length when they cross the declared origin.

Within a primary partition, FASTA and truth have the same deterministic record
order and a one-to-one ID relation. Similarity files are overlapping **views**:
a retained fragment can occur in `test`, `candidate_test`, and one stratum;
an excluded candidate occurs in `candidate_test` and one stratum but not strict
`test`. This intentional cross-file reuse is not train/test leakage.

## Root files

### `.chimera-bundle`

A small marker identifying a directory CHIMERA may replace when `--force` is
explicitly supplied. Its content is the bundle schema URN. Do not copy this
marker into an unrelated directory: it is a safety boundary, not decoration.

### `manifest.json`

The machine-facing bundle manifest contains:

| Key | Meaning |
|---|---|
| `schema` | `urn:chimera:benchmark-bundle:1`. |
| `tool` | CHIMERA name, exact version, `software_content_sha256`, `git_revision`, and nullable `git_dirty` provenance. The content receipt hashes the executable package sources and canonical schemas, so installed wheels remain identifiable when no Git checkout is present. |
| `data_model` | Alphabet; explicit linear/circular coordinate systems and semantics; opaque-header rule; source grouping; synthetic status. |
| `randomness` | Master seed, PRNG description, and semantic seed-derivation namespace. |
| `references` | Retained count, content-addressed input receipts, and preflight-exclusion count. Each `inputs` item has `content_id`, `role`, and `sha256`. |
| `splits` | Map from protocol name to the content also stored in each `split.json`. |

Input receipts contain no local path: `content_id` is `sha256:<digest>`, `role`
is `reference_fasta`, `reference_metadata`, or `external_similarity_table`, and
`sha256` is the same digest without its algorithm prefix.

### `resolved-config.json`

The fully validated semantic configuration after TOML and CLI overrides. Input
paths are replaced by sorted `sha256:<digest>` content IDs, `output_dir` is the
stable logical value `bundle`, split names are canonicalized, fragment lengths
are sorted, dates are serialized, and all defaults are materialized. The
operational `overwrite` flag is omitted because it does not define dataset
semantics.

### `schemas/`

An embedded snapshot of the JSON Schemas for the bundle, split manifest,
resolved configuration, and logical metadata/reference/sequence/truth/
assignment/exclusion rows. `chimera validate` requires the complete inventory
and verifies that it matches the schemas shipped by the validating CHIMERA
version. Use `chimera schema NAME` to print an installed schema.

### `execution.json`

Run-specific audit data: start/finish UTC timestamps, Python version, platform,
and completion status. It is intentionally excluded from `checksums.sha256` so
semantic files can reproduce while execution facts remain truthful.

### `references.tsv`

One row per retained source genome. Print the current header with
`chimera schema references`.

| Column | Meaning |
|---|---|
| `genome_id` | User grouping ID; joins assignments and truth `source_genome_id`. |
| `label` | `virus` or `host`. |
| `accession_version` | Common supplied accession version, otherwise empty. |
| `release_date` | Effective first-public-release date; grouped genomes use latest segment release. |
| `sequence_ids` | Canonical compact JSON array of member FASTA IDs in stable order. |
| `contig_count` | Number of member FASTA records. |
| `length_nt` | Sum of contig lengths. |
| `sha256` | `CHIMERA-GENOME-SHA256-v2` topology-aware whole-genome digest, invariant to contig order and strand orientation. |
| `source_input_ids` | Canonical compact JSON array of `sha256:<digest>` identifiers for the input FASTA files contributing member sequences. |
| `taxonomy` | Canonical compact JSON object of group-level `rank: value` pairs. |
| `metadata_extra` | Canonical compact JSON object; a source field appears here only when present with the same value on every grouped sequence. |

### `sequences.tsv` and `source-sequences.fasta.gz`

These files provide a lossless, one-row-per-retained-FASTA-record source
inventory. Their ordered `sequence_id` values must agree exactly. The source
FASTA contains normalized IUPAC DNA and identifiers only; `sequences.tsv`
contains the provenance required to authenticate fragment truth.

| `sequences.tsv` column | Meaning |
|---|---|
| `sequence_id` | Exact source FASTA identifier; joins truth `source_sequence_id`. |
| `genome_id` | Parent source genome; joins `references.tsv`. |
| `label` | `virus` or `host`. |
| `accession_version` | Sequence-level supplied accession/version, otherwise empty. |
| `release_date` | Exact sequence/segment first-public-release date, otherwise empty. |
| `topology` | Declared `linear` or `circular`. |
| `length_nt` | Normalized sequence length. |
| `sha256` | SHA-256 of the normalized forward sequence bytes. |
| `canonical_sha256` | `CHIMERA-CONTIG-SHA256-v2`: strand-invariant and topology-domain-tagged; circular records are also origin-rotation invariant. |
| `source_input_id` | `sha256:<digest>` identifier of the contributing input FASTA file. |
| `taxonomy` | Canonical compact JSON object of sequence-level taxonomy. |
| `metadata_extra` | Canonical compact JSON object preserving additional sequence-level metadata. |

### Root `excluded.tsv`

References removed during catalog preflight, currently same-class exact/content
duplicates accepted under `duplicate_policy = "drop"`. It uses the common
exclusion schema below with `split = reference_preflight` and fills
`duplicate_of`.

### `REPORT.md`

A human-readable run summary: tool/seed, reference counts, protocol test counts,
and interpretation limits. It is derived documentation, not a substitute for
the JSON manifests or TSV assignments.

### `checksums.sha256`

Sorted lowercase SHA-256 records in conventional form:

```text
<64 hexadecimal characters><two spaces><relative POSIX path>
```

Every regular bundle file is covered except `checksums.sha256` itself and
`execution.json`. Paths are relative and must remain within the bundle. Use
`chimera validate`, not an unchecked path-expanding parser, when accepting an
untrusted bundle.

## Protocol files

### `train.fasta.gz`, `test.fasta.gz`, and similarity views

Each record is exact synthetic DNA with an opaque ID of the form
`frag-` followed by a stable 32-character lowercase hexadecimal token. The ID
format is opaque: downstream code must not infer semantics from its spelling or
depend on its length. There is no description containing label, source, split,
date, or taxonomy.

`train` and `test` are the primary partitions. In 2A, every genome contributes
to both. In 2B–2E, source/content groups are disjoint.

For 2C:

- `candidate_test` is the entire genome-disjoint test proposal before the
  strict gate;
- strict `test` omits candidates whose best result passes both the configured
  above-identity and coverage gate;
- each candidate occurs in exactly one `test_strata` novelty view, irrespective
  of strict retention.

### `*.truth.tsv.gz`

One row per FASTA record. Print the current header with `chimera schema truth`.

| Column | Meaning |
|---|---|
| `sequence_id` | Opaque fragment ID; exact join key to FASTA. |
| `label` | Ground-truth `virus` or `host`. |
| `source_accession_version` | Exact source sequence/segment accession/version, possibly empty. |
| `source_genome_id` | Source `genome_id`; joins `references.tsv` and `assignments.tsv`. |
| `source_content_group_id` | Canonical source content group as `sha256:<whole-genome digest>`. |
| `source_sequence_id` | Original contig/segment FASTA ID. |
| `source_start` | Zero-based inclusive start on the forward source contig. |
| `source_end` | Zero-based exclusive unwrapped end on the forward source contig. |
| `coordinate_system` | `0-based-half-open` for linear sources or `0-based-half-open-circular` for circular sources. |
| `strand` | `+` if FASTA is the source interval, `-` if its reverse complement. |
| `fragment_length` | Exact nucleotide count; equals `source_end - source_start`. |
| `partition` | Semantic assignment of this observation: `train`, `test`, or `excluded`. |
| `view` | Serialized view: `train`, `test`, `candidate_test`, or `test_strata/<similarity_bin>`. |
| `similarity_bin` | 2C similarity bin; empty for other protocols/training. |
| `max_train_similarity` | Candidate's best train similarity on `[0,1]`; empty when unavailable/not applicable. |
| `nearest_train_genome_id` | Best matching training `genome_id`; empty when unavailable/not applicable. |
| `release_date` | Exact source sequence/segment release date, if supplied. |
| `synthetic` | Literal `true`; records are generated substrings. |

For a linear `strand = -` row, verify sequence truth as
`reverse_complement(source[source_start:source_end])`. For a circular row whose
end exceeds source length, first reconstruct the forward interval as
`source[source_start:] + source[:source_end - source_length]`, then apply the
strand. Coordinates never refer to the orientation of the emitted string.

`partition` is biological split semantics; `view` is file membership. A 2C
candidate rejected by the strict gate therefore has `partition = excluded` in
both `candidate_test` and its stratum view. Intentional duplication across
auxiliary views does not change that semantic assignment.

### `assignments.tsv`

One row per retained or protocol-excluded source genome. Test 2A records
`partition = both`; other protocols use `train`, `test`, or `excluded`.

| Column | Meaning |
|---|---|
| `genome_id` | Source genome join key. |
| `group_id` | Canonical content group (`sha256:<digest>`), or stable group ID. |
| `label` | `virus` or `host`. |
| `partition` | Final protocol disposition: `both`, `train`, `test`, or `excluded`. |
| `candidate_partition` | 2C pre-gate `train`/`test`; empty outside 2C. Excluded candidates retain `test`. |
| `reason` | Stable, auditable assignment/exclusion reason. |
| `release_date` | Source effective release date when present. |
| `taxon` | 2E value at the selected rank, otherwise empty. |
| `similarity_bin` | 2C candidate band, otherwise empty. |
| `nearest_train_genome_id` | 2C best training source, otherwise empty. |
| `max_train_similarity` | 2C best value on `[0,1]`, otherwise empty. |
| `similarity_coverage` | External aligned fraction if supplied; built-in screen leaves it empty. |
| `similarity_method` | External method string or built-in method descriptor. |
| `strict_gate_train_genome_id` | Training source for the hit that triggered strict exclusion; may differ from the numerical best hit. |
| `strict_gate_similarity` | Similarity of that exclusion-triggering hit. |
| `strict_gate_coverage` | Coverage of that hit; empty for the built-in screen. |
| `strict_gate_method` | Method descriptor for that hit. |

Use both `partition` and `candidate_partition` to distinguish a strict 2C
exclusion from a source never proposed for test.

### Per-protocol `excluded.tsv`

One row per source excluded from the final protocol. Empty files retain the
header.

| Column | Meaning |
|---|---|
| `genome_id` | Excluded source ID. |
| `label` | `virus` or `host`. |
| `split` | Canonical protocol name (`similarity`, `temporal`, etc.). |
| `reason` | Explicit exclusion reason. |
| `duplicate_of` | Retained representative for a root preflight duplicate; otherwise empty. |
| `source_sha256` | Canonical whole-genome source digest. |
| `source_accession_version` | Effective group accession/version, otherwise empty. |
| `release_date` | Effective source-genome release date, otherwise empty. |
| `nearest_train_genome_id` | Similarity exclusion's best training source, otherwise empty. |
| `max_train_similarity` | Similarity exclusion's best value, otherwise empty. |
| `similarity_coverage` | Coverage attached to that best value, otherwise empty. |
| `similarity_method` | Method attached to that best value, otherwise empty. |
| `strict_gate_train_genome_id` | Training source for the qualifying strict-gate hit. |
| `strict_gate_similarity` | Similarity of the qualifying strict-gate hit. |
| `strict_gate_coverage` | Coverage of the qualifying strict-gate hit, if applicable. |
| `strict_gate_method` | Method for the qualifying strict-gate hit. |

`missing_metadata = "exclude"` records temporal/taxonomic missingness here.
2C records candidates removed by its strict gate. Exclusions are part of the
scientific result and must be counted and archived.

### `split.json`

The stable split manifest contains:

| Key | Meaning |
|---|---|
| `schema` | `urn:chimera:split-manifest:1`. |
| `protocol` | Canonical name: `random`, `genome`, `similarity`, `temporal`, or `taxonomy`. |
| `protocol_id` | `2a` through `2e`. |
| `parameters` | Resolved protocol parameters, operators, grouping, cutoff/taxa, and method/source. |
| `validation` | Pass status and overlap/source counts; 2A is marked diagnostic-only. |
| `train`, `test` | Fragment statistics described below. |
| `truth_rows` | Train/test truth row counts. |
| `excluded_genomes` | Count written to per-split exclusion table. |
| `candidate_test` | 2C-only statistics for the complete candidate view. |

Fragment-statistic objects include `records`, `bases`, `gc_fraction`,
`ambiguous_fraction`, `records_by_label`, `records_by_length`,
`source_genomes`, and `records_by_genome`. Fractions are JSON `null` for an
empty view.

Protocol parameters are intentionally self-describing. For example, 2C records
the strict identity operator (`>`), strict coverage operator (`>=`), table or
built-in source, k, sketch size, thresholds, and similarity bands; 2D records its
inclusive cutoff and temporal semantics; 2E records rank and held-out taxa.

## Validation report counts

`chimera validate BUNDLE --json` returns the structured `ValidationReport`.
`primary_fasta_records_verified` and `primary_truth_rows_verified` count only
the train/test partitions summarized in `splits`. The corresponding
`auxiliary_*` fields count 2C `candidate_test` and all stratum views. Because a
candidate intentionally occurs in `candidate_test` and one stratum—and may also
occur in strict `test`—auxiliary counts are view-row counts, not unique fragment
or source counts.

`fasta_records_verified` and `truth_rows_verified` are compatibility totals,
each equal to primary plus auxiliary. The report also exposes
`assignment_rows_verified`, `checksums_verified`, the per-split primary
train/test counts, and the completed check names.

## Input metadata schema

The canonical discovery command prints:

```console
chimera schema metadata
```

Its columns are:

```text
sequence_id genome_id label accession_version release_date topology realm kingdom phylum class order family genus species
```

Only `sequence_id` is structurally required, but 2D and 2E require their
applicable values unless `missing_metadata = "exclude"`. See the
[user guide](USER_GUIDE.md#metadata) for aliases, grouping, and strict join
behavior.

## Compatibility policy

Bundle and split schema URNs identify the contract. Consumers should reject an
unknown major URN and must select fields by header name, not position. The v1
filesystem inventory is exact: extra entries are invalid even when checksum-
listed. Do not depend on row order as a randomization mechanism, opaque ID
internals, or prose in `REPORT.md`.
Resolve content IDs through a separately archived receipt/catalog when original
filenames matter. Record the CHIMERA version alongside the schema.
