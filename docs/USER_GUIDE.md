# CHIMERA user guide

This guide covers the command-line interface for CHIMERA 1.0.0. Run
`chimera COMMAND --help` for the installed version's authoritative flag list.

## Installation

CHIMERA supports Python 3.11–3.14. Install the current version from the GitHub
repository into an isolated environment:

```console
git clone https://github.com/Alexander-Mitrofanov/CHIMERA.git
cd CHIMERA
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
chimera --version
```

For development and contribution checks, install the development dependencies:

```console
python -m pip install -e '.[dev]'
make ci
```

See [CONTRIBUTING.md](../CONTRIBUTING.md) for the complete contributor setup.
The repository also provides a bounded `environment.yml` and a non-root
container recipe. The environment file is not a fully resolved lockfile;
archive an explicit environment lock/export or immutable container digest when
exact dependency reproduction is required.

## One command for Tests 2A–2E

The recommended interface is below. This path exists in a source checkout or
unpacked source distribution, not in an installed wheel:

```console
chimera suite --config examples/tiny/chimera.toml
```

`suite` always resolves all five evaluation protocols and rejects a narrowed
`splits` configuration or `--split` selection; use `generate` for a subset. It
performs reference discovery, normalization,
metadata joining, duplicate/conflict checks, split planning, fragment
generation, internal leakage checks, and atomic bundle publication. It does
not download references or contact a remote service.

To select protocols, use `generate` and repeat or comma-separate `--split`:

```console
chimera generate \
  --virus viruses.fna \
  --host hosts.fna \
  --metadata metadata.tsv \
  --outdir selected-benchmark \
  --split genome,similarity \
  --split taxonomy
```

Accepted names include `random`, `genome`, `similarity`, `temporal`,
`taxonomy`, and `all`; `2a` through `2e` and descriptive aliases are also
accepted.

## Preparing reference FASTA

Supply virus and host/non-viral references separately. Each option is
repeatable and accepts either a file or a directory:

```console
chimera inspect \
  --virus viral-a.fna --virus viral-directory/ \
  --host hosts.fna.gz \
  --metadata metadata.tsv
```

Directory discovery is recursive and deterministic. Supported suffixes,
matched case-insensitively, are `.fa`, `.fasta`, `.fna`, and their `.gz`
variants.

Input requirements:

- FASTA text must be UTF-8 and each record must have a non-empty sequence.
- The first whitespace-delimited header token is `sequence_id`. It must be
  globally unique across all virus and host files.
- A FASTA file containing more than one record requires `--metadata`, with one
  matching row per record. CHIMERA will not guess whether multiple records are
  independent genomes or segments/contigs of one genome. Without metadata,
  each accepted single-record FASTA file becomes one independent linear genome.
- IDs contain 1–255 ASCII characters, start with a letter or digit, and then
  use only letters, digits, `.`, `_`, `:`, `+`, or `-`.
- DNA is normalized to uppercase ungapped IUPAC symbols
  `ACGTRYSWKMBDHVN`. Whitespace is removed. `U`, gaps, and other characters are
  rejected instead of silently converted.
- Each genome must have at least one contig at least as long as every requested
  fragment length. Fragments never cross a contig boundary.
- A suite needs at least two independent virus content groups and two
  independent host content groups after integrity filtering, so every
  non-random split can retain both classes on both sides.

Source grouping uses the topology-aware `canonical_topology_aware_genome_sha256_v2`
digest; circular origin rotations are equivalent, while linear and circular
sources remain distinct. A separate topology-agnostic exact raw/RC audit makes
contradictory virus/host content fatal even when topology metadata differs.
Same-class duplicates within the same topology semantics either fail
(`duplicate_policy = "error"`, recommended while curating) or are removed
deterministically and recorded in root `excluded.tsv` (`duplicate_policy = "drop"`).

## Metadata

Metadata can be tab-separated UTF-8 or comma-separated when the filename ends
in `.csv`. Print the canonical columns with:

```console
chimera schema metadata
```

| Column | Requirement and meaning |
|---|---|
| `sequence_id` | **Required.** Exact FASTA ID; one row per FASTA record. |
| `genome_id` | Groups multiple contigs/segments as one indivisible source genome. Defaults to `sequence_id`. |
| `label` | Optional `virus` or `host`; if present it must agree with the input option used. |
| `accession_version` | Recommended stable accession including version, for example `NC_012345.2`. |
| `release_date` | ISO `YYYY-MM-DD` first public release date. Required by 2D unless missing records are explicitly excluded. |
| `topology` | `linear` (default) or `circular`, declared independently for every FASTA record. |
| `realm` … `species` | Canonical rank values used for provenance; the configured holdout rank is required for viral 2E records. |

The full recommended header is:

```text
sequence_id  genome_id  label  accession_version  release_date  topology  realm  kingdom  phylum  class  order  family  genus  species
```

The display above uses spaces for readability; the actual file must use tabs.
Rank aliases such as `tax_family` are accepted. Headers are stripped,
lowercased, and have spaces/hyphens converted to underscores. The ambiguous
legacy date headers `deposited_at` and `create_date` are rejected. Verify the
first-public date upstream and name the column exactly `release_date`.

Metadata joining is strict: every FASTA ID needs exactly one metadata row and
unused metadata rows are errors. Within one grouped genome, taxonomy at a rank
must be consistent. Its effective release date is the latest segment release,
because the complete grouped genome was not available earlier. A common
`accession_version` is retained only when all segments agree; otherwise the
group value is blank. Additional columns are accepted. Their per-record values
are preserved as JSON in `sequences.tsv` `metadata_extra`; a value is promoted
to `references.tsv` `metadata_extra` only when that key occurs with the same
value on every segment of the grouped genome.

`release_date` is not collection date. NCBI describes Virus `releaseDate` as
the first public release; record collection dates separately when they matter.
Pin `accession.version`, because sequence updates increment the version.

## TOML configuration

Configuration may use a `[benchmark]` table (recommended) or place its fields
at the TOML root. Relative paths in TOML are resolved from the configuration
file's directory. Relative paths passed as CLI flags are resolved from the
invocation directory. Explicit CLI flags override file values.

```toml
[benchmark]
schema_version = 1
virus_paths = ["references/viruses.fna.gz"]
host_paths = ["references/hosts.fna.gz"]
metadata_path = "references/metadata.tsv"
output_dir = "benchmark-v1"
splits = ["random", "genome", "similarity", "temporal", "taxonomy"]

seed = 42
test_fraction = 0.20
fragment_lengths = [150, 500]
fragments_per_genome = 200
strand_mode = "both"
max_ambiguous_fraction = 0.05
duplicate_policy = "error"
missing_metadata = "error"

temporal_cutoff = "2021-12-31"
taxonomy_rank = "family"
holdout_taxa = ["Exampleviridae"]
auto_holdout_count = 1

similarity_k = 21
sketch_size = 2000
max_train_similarity = 0.95
min_similarity_coverage = 0.85
similarity_bands = { high = 0.90, moderate = 0.70, low = 0.30 }
# similarity_table = "similarities.tsv"
```

Unknown keys and unsupported `schema_version` values fail. Hyphenated keys are
normalized to underscores, but the underscore form above is the stable style.

### Configuration reference

| TOML key | CLI flag | Default | Constraint/meaning |
|---|---|---:|---|
| `virus_paths` | `--virus` | required | File/directory list; repeat flag. |
| `host_paths` | `--host` | required | File/directory list; repeat flag. |
| `metadata_path` | `--metadata` | none | Needed for temporal/taxonomic suite protocols. |
| `output_dir` | `--outdir` | `benchmark-output` | Dedicated new bundle directory. |
| `splits` | `--split` | all | Protocol list. |
| `seed` | `--seed` | `42` | Integer master seed. |
| `test_fraction` | `--test-fraction` | `0.20` | Strictly between 0 and 1. |
| `fragment_lengths` | `--fragment-length` | `[500]` | Unique positive exact lengths; repeat flag. |
| `fragments_per_genome` | `--fragments-per-genome` | `100` | At least twice the number of requested lengths, so every genome/length stratum can enter both 2A partitions. |
| `strand_mode` | `--strand` | `both` | `both` or `forward`. |
| `max_ambiguous_fraction` | `--max-ambiguous-fraction` | `0.05` | Allowed non-ACGT fraction in each emitted fragment. |
| `duplicate_policy` | `--duplicate-policy` | `error` | `error` or audited `drop`. |
| `missing_metadata` | `--missing-metadata` | `error` | `error` or explicit per-protocol `exclude`. |
| `temporal_cutoff` | `--release-date-cutoff` | auto | Inclusive ISO date; auto mode selects a class-viable cutoff near the target fraction. |
| `taxonomy_rank` | `--holdout-rank` | `family` | One supplied metadata rank name. |
| `holdout_taxa` | `--holdout-taxon` | stable auto | Explicit viral value(s); repeat flag. |
| `auto_holdout_count` | `--auto-holdout-count` | `1` | Positive count when taxa are auto-selected. |
| `similarity_table` | `--similarity-table` | built-in | Versioned external all-candidate-vs-train TSV. |
| `similarity_k` | `--similarity-k` | `21` | Odd canonical k-mer length from 5 through 63. |
| `sketch_size` | `--sketch-size` | `2000` | Bottom-k size, at least 100. |
| `max_train_similarity` | `--max-train-similarity` | `0.95` | Strict identity maximum; must be at least the high-band boundary. |
| `min_similarity_coverage` | `--min-similarity-coverage` | `0.85` | Minimum external aligned fraction for an above-threshold hit to exclude. |
| `similarity_bands` | three `--similarity-*` flags | `.90/.70/.30` | Must satisfy `1 >= high > moderate > low >= 0`. |
| `overwrite` | `--force` | `false` | Replace only a recognized CHIMERA bundle via atomic backup/commit. |

The CLI band flags are `--similarity-high`, `--similarity-moderate`, and
`--similarity-low`.

## Fragment generation

For each genome, CHIMERA balances the requested lengths as evenly as possible
within the configured record count. Coordinates are sampled with replacement;
all eligible start positions across eligible contigs have equal probability.
Consequently duplicate fragment sequence or coordinates are possible by
design, but fragment IDs remain unique. No fragment spans contigs.

With `strand_mode = "both"`, each accepted coordinate is deterministically
assigned forward (`+`) or reverse-complement (`-`) orientation. `forward`
disables reverse strands. Truth coordinates always refer to the forward source
contig. Linear records use `[source_start, source_end)` with `source_end` no
greater than source length. Circular records use an unwrapped forward interval:
`source_start` is within the source, `source_end = source_start +
fragment_length`, and `source_end` may exceed source length when the fragment
wraps the declared origin.

Coordinates exceeding `max_ambiguous_fraction` are rejected with a bounded
attempt budget. If a genome cannot provide the requested accepted records,
CHIMERA stops with an actionable error; it never silently emits a short
dataset. This is substring sampling, not an empirical read model.

## Understanding each protocol

### Test 2A — random fragment diagnostic

CHIMERA first generates fragments for every source genome, then performs a
deterministic shuffled split independently within each genome. Each genome has
non-empty train and test contributions, fragment IDs never overlap, and source
genome overlap is expected. The combined output order is deterministically
shuffled across classes and sources.

Use 2A to detect basic pipeline failures and compare with older random-split
work. Do not use it as the headline result for unseen-genome discovery.

### Test 2B — genome holdout

Whole genomes are assigned before fragments are generated. Assignment is
label-stratified toward `test_fraction`. Canonical whole-genome SHA-256 groups
are invariant to contig input order and reverse-complement representation, so
equivalent content cannot cross partitions. This tests unseen source genomes,
not necessarily unseen homologous families.

### Test 2C — similarity-filtered holdout

CHIMERA makes a label-stratified, genome/content-disjoint candidate split and
compares every candidate test genome with training. It writes two primary test
views:

- `candidate_test` is the complete proposal, before the strict gate.
- `test` is the strict retained set.

For an external-table result, an above-threshold candidate is excluded only
when `similarity > max_train_similarity` **and**
`coverage >= min_similarity_coverage`. Equality at the identity threshold is
retained. A reported above-threshold hit below minimum coverage stays in the
strict set with an explicit reason. Detected external hits always include
aligned coverage. Built-in MinHash has no alignment coverage, so its coverage
is absent and the identity gate applies.

Novelty strata are computed for the complete candidate set:

- `high_similarity`: similarity `>= high`;
- `moderate_similarity`: `>= moderate` and `< high`;
- `low_similarity`: `>= low` and `< moderate`;
- `distant_detectable`: a measured similarity `< low`;
- `no_detectable_match`: no shared canonical k-mer or external detection.

Every stratum is written even when empty. The strata do not redefine strict
membership. Report strict performance and candidate novelty curves separately.

#### Built-in screen

The offline engine hashes canonical unambiguous k-mers, never joins k-mers
across contig boundaries, retains a stable bottom-k sketch, estimates Jaccard,
and converts it to a Mash-style identity estimate. It is exact only when the
sketch retains every unique canonical k-mer. Ambiguity-bearing k-mers are
ignored. No shared k-mer is represented as no detectable match.

This value is not alignment-derived ANI and is not a validated universal viral
taxonomic threshold. A small or poorly chosen k, limited sketch, fragmented
assembly, repeats, and divergent genome length can alter the estimate. For a
publication's decisive similarity analysis, use a biologically appropriate,
versioned alignment workflow and archive its complete input table.

#### External similarity table

Pass a UTF-8 tab-separated file:

```text
query_genome_id\treference_genome_id\tsimilarity\tcoverage\tcoverage_definition\tmethod
candidate-1\ttrain-7\t0.9412\t0.9120\taligned_fraction_shorter\tskani-0.x
candidate-1\ttrain-9\t0.8120\t0.7340\taligned_fraction_shorter\tskani-0.x
candidate-2\ttrain-7\t\t\taligned_fraction_shorter\tskani-0.x
candidate-2\ttrain-9\t\t\taligned_fraction_shorter\tskani-0.x
```

All six columns are required. Similarity and coverage are fractions on `[0,1]`,
not percentages. References must be CHIMERA training `genome_id` values and
queries must be candidate test IDs. The table must contain exactly one row for
every candidate-by-training pair. For a nondetected pair, keep both IDs and
leave both `similarity` and `coverage` blank; leaving only one blank is invalid.
`coverage_definition` must be `aligned_fraction_shorter`, and `method` must
contain a tool/version identifier, including for nondetected pairs. Mixed
queries select the greatest detected similarity; a query whose complete row
set is nondetected becomes `no_detectable_match` and has no arbitrary nearest
training ID. Similarity ties are broken by lexical reference ID.

Because candidate assignment depends on the configured seed and inputs, first
generate a preliminary built-in-similarity bundle and read
`2c_similarity_filtered/assignments.tsv` for rows whose
`candidate_partition` is `train` or `test`. Compute the external all-pairs table
for those IDs, then generate a **new** final output directory with the same
references, seed, and split parameters plus `similarity_table`. The proposal is
deterministic and independent of the similarity source. `--dry-run` confirms
plan counts but does not write the assignment IDs. Archive both the external
table and its tool/version, database snapshot, command, identity definition,
and coverage definition.

### Test 2D — temporal holdout

With an explicit cutoff, training contains genomes with effective
`release_date <= cutoff`; test contains genomes with a later date. The cutoff
is inclusive. Missing dates either fail or appear in per-split exclusions,
never in training, according to `missing_metadata`.

If no cutoff is supplied, CHIMERA deterministically selects a date that leaves
at least one virus and one host on both sides and approximates the target test
fraction. For sequence-content-equivalent groups, the earliest represented
public availability is used to prevent a later duplicate representation from
crossing the cutoff. For a multi-contig genome grouped by one `genome_id`, its
metadata date is the latest segment release.

A run against today's sequences and taxonomy is a release-date-filtered
**retrospective** split. It does not reconstruct what a researcher knew at the
cutoff. A prospective claim requires an archived historical database and
taxonomy snapshot, pinned accession versions, and documented preprocessing.
CHIMERA's CLI does not manufacture or certify such a snapshot.

### Test 2E — taxonomic holdout

All viral genomes with selected values at `taxonomy_rank` are assigned to test;
other represented viral taxa train the classifier. Hosts receive an independent
label-stratified genome/content-disjoint split so both classes occur on both
sides. Explicit taxa are recommended for a preregistered analysis; otherwise a
stable seeded rule selects `auto_holdout_count` taxa while preserving at least
one viral training taxon.

Matching is case-insensitive but otherwise exact on the supplied string.
CHIMERA does not resolve taxonomic identifiers, ancestors, spelling variants,
synonyms, or renamed taxa. Normalize against a versioned taxonomy upstream.
The data must contain at least two represented viral taxa at the chosen rank,
and the holdout cannot consume all of them. Missing rank values fail or are
explicitly excluded. A held-out family can still be highly similar to training,
so do not treat 2E as a substitute for 2C.

## Safe, deterministic output

Generation occurs in a sibling staging directory. Files are flushed and
atomically replaced, then the complete staged directory is committed. CHIMERA
refuses `/`, the current directory, and the user's home directory as an output
target. An existing directory is never overwritten without `--force`, and an
unrecognized directory is never replaced even with `--force`.

Random decisions use Python `random.Random` with domain-separated,
semantic BLAKE2b-derived sub-seeds. Input ordering does not alter semantic
assignment. Tables are stably ordered, JSON is sorted and rejects non-finite
numbers, and gzip has no filename and timestamp zero. `checksums.sha256` covers
semantic bundle files except itself and `execution.json`; the latter contains
run timestamps and platform details. Treat exact repeatability as conditional
on identical bytes, resolved configuration, CHIMERA version, and compatible
runtime behavior; always archive those facts.

## Inspection, validation, and schemas

Preflight and print a normalized source inventory without generating:

```console
chimera inspect --virus viruses.fna --host hosts.fna --metadata metadata.tsv
chimera inspect --virus viruses.fna --host hosts.fna --metadata metadata.tsv --json
```

`--duplicate-policy drop` is also available for inspection. The human table
prints `genome_id`, `label`, length, contig count, release date, and canonical
SHA-256.

Resolve the complete suite without output:

```console
chimera suite --config benchmark.toml --dry-run
```

Validate a completed bundle independently:

```console
chimera validate benchmark-v1
chimera validate benchmark-v1 --json
```

Discover stable public table columns:

```console
chimera schema metadata
chimera schema truth
chimera schema references
chimera schema sequence-row
```

Every bundle also carries `source-sequences.fasta.gz` and a one-row-per-record
`sequences.tsv` inventory. Their ordered IDs must agree; the inventory records
the source `genome_id`, topology, exact and canonical hashes, sequence-level
date/taxonomy, content-addressed `source_input_id`, and extra metadata used to
authenticate truth. `references.tsv` aggregates those receipts as the JSON
`source_input_ids` array. Bundle manifests and resolved configuration likewise
use `sha256:<digest>` content IDs rather than disclosing absolute local paths;
archive a separate controlled mapping when original filenames matter.

`chimera validate --json` reports primary FASTA/truth counts for the train/test
partitions separately from auxiliary counts for the overlapping 2C
`candidate_test` and stratum views. The aggregate
`fasta_records_verified`/`truth_rows_verified` values are sums of primary and
auxiliary counts, so they are not counts of unique biological observations.

Exit codes are `0` for success, `2` for configuration/input/user errors, and
`3` for integrity-validation failure. Global options precede the subcommand;
for example, `chimera --quiet suite --config benchmark.toml`. `--quiet`
suppresses progress messages, while
`--log-level DEBUG|INFO|WARNING|ERROR` controls diagnostics on standard error.

## Migrating from the legacy executable

The old console name remains a temporary compatibility shim:

```console
metagenome-generator suite --config benchmark.toml
```

It prints `WARNING: 'metagenome-generator' is deprecated; use 'chimera'
instead.` and delegates to the same parser. Replace the executable name now.

The rewritten tool intentionally does not guess legacy option spellings or
promise the historical flat output layout. A safe migration is:

1. Build a `schema_version = 1` TOML file with explicit input, output, seed,
   lengths, split, temporal, taxonomy, and similarity choices.
2. Add a one-to-one metadata table, including pinned `accession_version`,
   release dates, and version-normalized taxonomy.
3. Run `chimera inspect`, then `chimera suite --dry-run`.
4. Generate into a new directory; do not point `--force` at legacy data.
5. Update consumers to join opaque FASTA IDs with `*.truth.tsv.gz`, and obtain
   source assignment from `assignments.tsv`, rather than parsing headers.
6. Compare record/class/source counts and register the new bundle as a new
   dataset version; it is not byte- or split-compatible with legacy output.

## Common failures

- **“requires at least two independent … content groups”**: add genuinely
  distinct sources for that class; duplicated accessions do not create a valid
  holdout.
- **Genome cannot emit a requested length**: shorten the fragment length or
  replace incomplete references. CHIMERA never crosses contigs.
- **Bounded ambiguous-fragment rejection**: improve/reference-mask policy,
  reduce length, or explicitly raise `max_ambiguous_fraction`; do not silently
  accept fewer records.
- **Temporal split not viable**: add date-diverse records for both labels or
  choose an explicit cutoff with both classes on each side.
- **Taxonomy split not viable**: supply at least two normalized viral values at
  the rank and leave one in training.
- **Strict similarity loses a class**: lower candidate relatedness, adjust a
  scientifically justified gate, or add independent sources. Do not accept an
  invalid class-empty benchmark.
- **Output already exists**: choose a new versioned directory. Use `--force`
  only for a recognized disposable/rerunnable CHIMERA bundle.

See [output formats](OUTPUT_FORMATS.md) before writing downstream parsers,
review the [methodology](METHODOLOGY.md), and complete the
[dataset datasheet](../DATASHEET.md) with dataset-specific evidence.
