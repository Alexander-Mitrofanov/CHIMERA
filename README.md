# CHIMERA

CHIMERA creates deterministic, leakage-aware virus-versus-host sequence
classification benchmarks. One command produces the basic random-fragment
diagnostic and four complementary generalization tests, with fragment-level truth,
source assignments, exclusions, checksums, and resolved provenance.

> **Scope.** CHIMERA samples exact DNA substrings from reference genomes. It
> does not simulate platform errors, abundance profiles, paired-end libraries,
> or complete metagenomic communities. The generated records are synthetic
> classification fragments, not empirical sequencing reads.

## Quick start

CHIMERA supports Python 3.11–3.14. Install the current version from this
repository into an isolated environment:

```console
git clone https://github.com/Alexander-Mitrofanov/CHIMERA.git
cd CHIMERA
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
chimera --version
chimera --help
```

Run the bundled tiny fixture without network access:

```console
chimera suite --config examples/tiny/chimera.toml
```

That exact command creates all Tests 2A–2E in
`examples/tiny/tiny-benchmark/`. The configuration contains two synthetic
viruses in two fictional families and two synthetic hosts, with release dates
on both sides of the cutoff. See [`examples/tiny/README.md`](examples/tiny/README.md).

For your own references, start with:

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

`suite` always selects all five protocols; use `generate --split ...` for a
subset. Explicit command-line flags override TOML values. Use `--dry-run` to
load, normalize, integrity-check, and resolve every split without writing a
bundle.

## Evaluation suite

| Test | Stable directory | Question and assignment unit |
|---|---|---|
| **2A — Random fragment** | `2a_random_fragment/` | Basic diagnostic. Fragments from every genome are shuffled within genome into train and test, so genome overlap is deliberate. |
| **2B — Genome holdout** | `2b_genome_holdout/` | Can a model recognize unseen genomes? Whole topology-aware, reverse-complement-invariant digest-v2 groups are assigned before fragment generation. |
| **2C — Similarity filtered** | `2c_similarity_filtered/` | Does performance degrade as references become less similar to training? A genome-disjoint candidate test is scored against training, stratified, and subjected to a strict gate. |
| **2D — Temporal holdout** | `2d_temporal_holdout/` | How does a model perform on accessions first released after a cutoff? Training uses `release_date <= cutoff`; test uses later dates. |
| **2E — Taxonomic holdout** | `2e_taxonomic_holdout/` | Can a model generalize to viral taxa absent from training? Selected viral rank values are held out in full. |

Test 2A is useful for continuity with conventional random evaluation, but it is
not the decisive generalization result: homologous biological entities can
leak across random record splits. Report 2B–2E separately, including exclusions
and per-similarity-stratum results. See the [methodology](docs/METHODOLOGY.md).

### Strict and candidate similarity sets

Test 2C deliberately writes both views:

- `candidate_test.*` contains every genome in the original genome-disjoint test
  proposal, including candidates rejected by the strict similarity gate.
- `test.*` is the strict primary set. A candidate is excluded when **any**
  reported external training hit is **greater than** `max_train_similarity`
  with coverage at least `min_similarity_coverage`; identity equality is
  retained. Built-in MinHash has no alignment coverage, so its identity gate
  applies directly.
- `test_strata/*` partitions the complete candidate set into `high_similarity`,
  `moderate_similarity`, `low_similarity`, `distant_detectable`, and
  `no_detectable_match`. These strata therefore remain
  useful for a performance-versus-novelty curve even when a record is absent
  from strict `test.*`.

The dependency-free screen is a canonical k-mer bottom-k MinHash calculation
reported as a Mash-style identity estimate. It is not alignment-derived ANI,
and its output must not be interpreted as a universal viral taxonomic boundary.
For publication analyses, supply a versioned all-candidate-versus-train table
from an appropriate validated alignment workflow with `--similarity-table`.

### Temporal and taxonomy qualifications

`release_date` means the accession's first public release date, not its sample
collection date. A normal CHIMERA 2D run uses current reference content and
taxonomy, so it is a **release-date-filtered retrospective** evaluation. Call it
prospective only if the entire reference and taxonomy inputs are a documented,
immutable historical snapshot available at the cutoff.

Taxonomic holdout compares the supplied rank strings case-insensitively. It
does not resolve synonyms, lineages, merged taxon identifiers, or spelling
variants. Normalize taxonomy upstream and retain the versioned source table.
Taxonomic novelty and sequence novelty are different claims; report 2C and 2E
independently.

## Auditable outputs

Each bundle includes:

```text
benchmark-v1/
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
├── 2b_genome_holdout/
├── 2c_similarity_filtered/
├── 2d_temporal_holdout/
└── 2e_taxonomic_holdout/
```

Every protocol has compressed train/test FASTA and truth tables,
`assignments.tsv`, `excluded.tsv`, and `split.json`. FASTA identifiers are
opaque and contain no class or source information. `sequences.tsv` and
`source-sequences.fasta.gz` form a one-to-one, per-source-sequence inventory for
truth reconstruction. Linear truth coordinates are conventional 0-based,
half-open forward-source intervals. Circular intervals use an unwrapped
0-based, half-open representation whose `source_end` may exceed source length
when the fragment crosses the declared origin; reverse-strand coordinates still
refer to the forward source. The complete tree and schemas are in
[`docs/OUTPUT_FORMATS.md`](docs/OUTPUT_FORMATS.md).

Generation is staged beside the destination and committed atomically. Existing
output is never replaced unless `--force` is given, and even then CHIMERA only
replaces a directory carrying its bundle marker. Semantic files use stable
ordering, canonical JSON/TSV, gzip timestamp zero, and BLAKE2-derived semantic
sub-seeds. `execution.json` contains run-time facts and is intentionally omitted
from `checksums.sha256`; with identical resolved inputs, configuration, CHIMERA
version, and compatible Python behavior, checksummed artifacts are repeatable.
Input receipts and resolved configuration use `sha256:<digest>` content IDs,
not absolute local paths.

## Inspect, validate, and discover schemas

```console
chimera inspect --virus viruses.fna --host hosts.fna --metadata metadata.tsv
chimera inspect --virus viruses.fna --host hosts.fna --metadata metadata.tsv --json
chimera validate benchmark-v1
chimera validate benchmark-v1 --json
chimera schema metadata
chimera schema truth
chimera schema references
chimera schema sequence-row
```

The JSON validation report separates primary train/test records from auxiliary
2C candidate/stratum views, which intentionally repeat records across views.

Exit status `0` means success, `2` means a configuration/input/user error, and
`3` means bundle-integrity validation failed. Diagnostics go to standard error;
machine-facing generation results go to standard output as JSON.

## Documentation

- [User guide](docs/USER_GUIDE.md): installation, FASTA/metadata preparation,
  flags, TOML, external similarities, and troubleshooting.
- [Methodology](docs/METHODOLOGY.md): algorithms, assumptions, scientific
  interpretation, and references.
- [Output formats](docs/OUTPUT_FORMATS.md): exact bundle tree and columns.
- [Dataset datasheet](DATASHEET.md): intended uses, risks, and maintenance.
- [Contributing](CONTRIBUTING.md), [security policy](SECURITY.md), and
  [changelog](CHANGELOG.md).

## Migration from MetagenomeGenerator

The maintained executable is `chimera`. The compatibility executable
`metagenome-generator` invokes the same CLI but emits a deprecation warning.
Migrate scripts now:

```diff
- metagenome-generator suite --config benchmark.toml
+ chimera suite --config benchmark.toml
```

Legacy ad-hoc output layouts and option names are not silently inferred. Create
a schema-versioned TOML file, inspect normalized inputs, run `--dry-run`, and
consume the manifest/truth tables rather than parsing FASTA headers. Detailed
mapping guidance is in the [user guide](docs/USER_GUIDE.md#migrating-from-the-legacy-executable).

## Citation, license, and support

CHIMERA 1.0.0 does not claim a project DOI. Cite the exact software version and
repository URL using [`CITATION.cff`](CITATION.cff), and cite each deposited
benchmark dataset with its own persistent identifier. Do not reuse a methods
paper DOI as a software or dataset DOI.

CHIMERA is released under the [MIT License](LICENSE). Report vulnerabilities
through the private process in [SECURITY.md](SECURITY.md); use the repository
issue tracker for reproducible non-security defects.
