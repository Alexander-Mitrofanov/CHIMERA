# CHIMERA

CHIMERA builds reproducible benchmark datasets for DNA classifiers that
distinguish viral sequence from host or other non-viral sequence. It takes
labeled reference genomes, samples exact fixed-length fragments, and writes
train/test partitions with machine-readable truth, source provenance,
exclusions, configuration, schemas, and checksums.

CHIMERA is a **dataset generator and validator**. It does not train a model,
predict labels, download reference databases, or simulate sequencing-platform
errors. It is also not a whole-community metagenome simulator: its generated
records are controlled DNA fragments sampled from the supplied references.

## Why use CHIMERA?

A random split of fragments can place sequence from the same source genome in
both training and test data. That setup is useful as a pipeline diagnostic, but
it can greatly overstate performance on novel biological material. CHIMERA
therefore provides five complementary partitioning protocols:

| Protocol | `--split` value | What it measures |
|---|---|---|
| Random fragment | `random` | A basic diagnostic in which every source genome contributes fragments to both partitions. It does **not** measure unseen-genome generalization. |
| Genome holdout | `genome` | Performance on source genomes that are absent from training. Whole content-equivalent genome groups are assigned before fragments are generated. |
| Similarity filtered | `similarity` | Performance as test genomes become less similar to training genomes. Candidate genomes are stratified by similarity and can be removed by a strict identity/coverage gate. |
| Temporal holdout | `temporal` | Retrospective performance on accessions first released after an inclusive cutoff date. |
| Taxonomic holdout | `taxonomy` | Generalization to selected viral taxa whose supplied rank values are absent from training. |

`chimera suite` generates all five protocols in one bundle. Use
`chimera generate --split ...` when only selected protocols are needed.

## Installation

CHIMERA supports Python 3.11–3.14. Install the current version from this
repository in an isolated environment:

```console
git clone https://github.com/Alexander-Mitrofanov/CHIMERA.git
cd CHIMERA
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
chimera --version
```

## Try the bundled example

The repository includes a small, entirely synthetic dataset that requires no
network access:

```console
chimera suite --config examples/tiny/chimera.toml
chimera validate examples/tiny/tiny-benchmark
```

The example contains two fictional viral genomes and two fictional host
genomes. It creates a complete bundle under
`examples/tiny/tiny-benchmark/`. See
[`examples/tiny/README.md`](examples/tiny/README.md) for details.

## Use CHIMERA with your own references

### 1. Prepare viral and host FASTA

Supply viral references with `--virus` and host/non-viral references with
`--host`. Each option accepts a FASTA file or a directory and may be repeated.
CHIMERA discovers `.fa`, `.fasta`, and `.fna` files, including `.gz` variants.
Sequences must use IUPAC DNA symbols.

Inputs are treated as untrusted data. CHIMERA validates FASTA structure,
normalizes sequence case, rejects conflicting identifiers and labels, and
checks for exact cross-class content conflicts before generating fragments.

### 2. Provide sequence metadata

A UTF-8 tab-separated metadata file is strongly recommended and is required
for multi-record FASTA inputs, temporal holdout, and taxonomic holdout. The
`sequence_id` column must match the FASTA record identifier. Use `genome_id` to
group contigs or genome segments that must never be split independently.

Typical columns are:

```text
sequence_id	genome_id	accession_version	release_date	family	genus	topology
virus_contig_1	virus_genome_1	V_000001.1	2019-04-01	Exampleviridae	Examplevirus	linear
host_contig_1	host_genome_1	H_000001.1	2018-07-10			linear
```

- `release_date` is the accession's first public release date in ISO
  `YYYY-MM-DD` form; it is not a collection date.
- Taxonomy values are supplied by the user. CHIMERA compares them
  case-insensitively but does not resolve names, lineages, synonyms, or taxon
  identifiers.
- `topology` is `linear` or `circular` and controls eligible fragment
  coordinates and content hashing.
- `label` may be included, but it must agree with the `--virus` or `--host`
  input channel.

Print the complete metadata header contract with:

```console
chimera schema metadata
```

### 3. Inspect references before generation

```console
chimera inspect \
  --virus references/viruses.fna \
  --host references/hosts.fna \
  --metadata references/metadata.tsv
```

Add `--json` for a machine-readable inventory. Inspection reports normalized
genome IDs, labels, lengths, contig counts, effective release dates, and
canonical content hashes without generating fragments.

### 4. Generate the complete benchmark suite

```console
chimera suite \
  --virus references/viruses.fna \
  --host references/hosts.fna \
  --metadata references/metadata.tsv \
  --outdir benchmark-v1 \
  --fragment-length 150 \
  --fragment-length 500 \
  --fragments-per-genome 200 \
  --seed 42
```

Important behavior:

- `--fragment-length` is repeatable. CHIMERA balances each genome's requested
  fragment count across the configured lengths.
- Coordinates are sampled uniformly with replacement from eligible positions;
  duplicate sequences or coordinates are therefore possible by design.
- `--strand both` samples forward and reverse-complement orientations.
- `--dry-run` validates inputs and resolves every split without writing a
  bundle.
- CHIMERA refuses to replace an existing output directory unless `--force` is
  supplied, and even then it only replaces a recognized, valid CHIMERA bundle.

### 5. Validate the completed bundle

```console
chimera validate benchmark-v1
chimera validate benchmark-v1 --json
```

Validation independently checks the bundle layout, schemas, checksums, source
inventory, reconstructed fragment sequences and coordinates, assignment
completeness, counts, exclusions, and protocol-specific leakage invariants.
Exit status `0` means success, `2` means an input/configuration error, and `3`
means bundle integrity validation failed.

## Generate only selected protocols

Use `generate` and repeat or comma-separate `--split` values:

```console
chimera generate \
  --virus references/viruses.fna \
  --host references/hosts.fna \
  --metadata references/metadata.tsv \
  --outdir genome-and-similarity \
  --split genome,similarity \
  --fragment-length 500 \
  --fragments-per-genome 100 \
  --seed 42
```

Accepted values are `random`, `genome`, `similarity`, `temporal`, `taxonomy`,
and `all`.

## Configuration files

Every generation option can be stored in TOML. Paths are resolved relative to
the configuration file:

```toml
[benchmark]
schema_version = 1
virus_paths = ["references/viruses.fna"]
host_paths = ["references/hosts.fna"]
metadata_path = "references/metadata.tsv"
output_dir = "benchmark-v1"
splits = ["random", "genome", "similarity", "temporal", "taxonomy"]

seed = 42
test_fraction = 0.20
fragment_lengths = [150, 500]
fragments_per_genome = 200
strand_mode = "both"

temporal_cutoff = "2021-12-31"
taxonomy_rank = "family"
holdout_taxa = ["Exampleviridae"]
```

Run it with:

```console
chimera suite --config benchmark.toml
```

Explicit command-line flags override TOML values. The
[user guide](docs/USER_GUIDE.md) documents every option, constraint, and
default.

## Bundle contents

A complete suite has the following top-level structure:

```text
benchmark-v1/
├── .chimera-bundle
├── manifest.json
├── resolved-config.json
├── execution.json
├── references.tsv
├── sequences.tsv
├── source-sequences.fasta.gz
├── excluded.tsv
├── REPORT.md
├── checksums.sha256
├── schemas/
├── random_fragment/
├── genome_holdout/
├── similarity_filtered/
├── temporal_holdout/
└── taxonomic_holdout/
```

Each protocol directory contains compressed train/test FASTA, matching truth
tables, source assignments, exclusions, statistics, and a split manifest. FASTA
identifiers are opaque: labels and source identities are available only through
the truth tables, preventing accidental label leakage through record names.

The similarity-filtered protocol also writes the complete pre-gate candidate
view and similarity strata. Its built-in canonical k-mer MinHash calculation
is an offline screening estimate, not alignment-derived ANI. For analyses that
require alignment identity and coverage, provide a complete, versioned
candidate-versus-training table with `--similarity-table`.

See [output formats](docs/OUTPUT_FORMATS.md) for the exact directory contract
and column definitions.

## Scientific interpretation

- Treat random-fragment results as a diagnostic, not evidence of
  unseen-genome recognition.
- Genome holdout prevents shared source genomes and equivalent whole-genome
  content across partitions, but it does not guarantee low homology.
- Similarity values and thresholds are study choices, not universal viral
  species boundaries.
- A temporal run over present-day references is retrospective. A prospective
  claim requires an immutable historical sequence and taxonomy snapshot.
- Taxonomic holdout depends on the supplied normalized rank strings and is not
  a substitute for sequence-similarity analysis.
- CHIMERA validates mechanics and declared invariants; it cannot certify the
  biological correctness of user-provided labels, dates, taxonomy, or external
  similarity evidence.

The detailed algorithms and assumptions are documented in
[Methodology](docs/METHODOLOGY.md).

## Reproducibility and safety

CHIMERA derives random decisions from a master seed and stable semantic
identities rather than input order. It writes canonical JSON/TSV, deterministic
gzip streams, content-addressed input receipts, the fully resolved
configuration, software provenance, and SHA-256 checksums. With identical input
bytes, configuration, CHIMERA version, and compatible Python behavior,
checksummed outputs are repeatable.

Generation occurs in a staging directory and is committed atomically. Absolute
local input paths are not stored in checked provenance, and broad or
unrecognized output directories are never deleted.

## Documentation and support

- [User guide](docs/USER_GUIDE.md): complete CLI/configuration reference and
  troubleshooting.
- [Methodology](docs/METHODOLOGY.md): algorithms, leakage controls,
  assumptions, and scientific limitations.
- [Output formats](docs/OUTPUT_FORMATS.md): bundle layout, table columns, and
  schemas.
- [Dataset datasheet](DATASHEET.md): template for documenting a generated
  benchmark dataset.
- [Contributing](CONTRIBUTING.md), [security policy](SECURITY.md), and
  [changelog](CHANGELOG.md).

The maintained executable is `chimera`. The legacy
`metagenome-generator` executable delegates to the same CLI with a deprecation
warning.

CHIMERA is released under the [MIT License](LICENSE). Cite the exact software
version and repository URL using [`CITATION.cff`](CITATION.cff), and cite each
deposited benchmark dataset with its own persistent identifier. Report security
issues privately through [SECURITY.md](SECURITY.md); use the GitHub issue
tracker for reproducible non-security defects.
