# Methodology and scientific interpretation

## Purpose

CHIMERA constructs controlled binary sequence datasets for evaluating whether a
classifier distinguishes virus-derived fragments from host/non-viral fragments
under increasingly demanding generalization protocols. It provides source
truth and leakage checks; it does not train or score the classifier.

The central design principle is that the unit sampled into train or test should
match the biological claim. A random fragment result asks whether a model can
recognize new fragments of already seen genomes. It does not answer whether the
model recognizes unseen genomes, distant sequence, later releases, or excluded
taxa. Random record splits can inflate estimates when homologous entities occur
on both sides, as discussed in a 2024 Nature Methods analysis
([DOI 10.1038/s41592-024-02362-y](https://doi.org/10.1038/s41592-024-02362-y)).
Leakage-aware partitioning is also the focus of DataSAIL
([DOI 10.1038/s41467-025-58606-8](https://doi.org/10.1038/s41467-025-58606-8))
and sequence-identity-aware SpanSeq
([DOI 10.1093/nargab/lqae106](https://doi.org/10.1093/nargab/lqae106)).
These works motivate CHIMERA's study design; their publication does not by
itself validate CHIMERA or make its protocols interchangeable with those tools.

## Reference data model and integrity boundary

One FASTA record is a contig or segment identified by `sequence_id`. Metadata
may map several records to one `genome_id`; that genome is the indivisible
biological source for Tests 2B–2E. The configured input channel supplies the
binary `virus` or `host` label, and an optional metadata label must agree. A
multi-record FASTA always requires an exhaustive metadata table because record
boundaries alone do not say whether records are independent genomes or grouped
segments. Without metadata, each accepted single-record FASTA file is treated
as an independent linear genome.

Sequences are normalized to uppercase ungapped IUPAC DNA. Each contig receives
a `CHIMERA-CONTIG-SHA256-v2` digest that is strand-invariant and explicitly
domain-tagged as linear or circular. Circular contigs are additionally invariant
to rotation of the declared origin. The `CHIMERA-GENOME-SHA256-v2` digest treats
these contig digests as a multiset, so order and names do not create a distinct
source while duplicate contigs remain significant. A separate
`CHIMERA-TOPOLOGY-AGNOSTIC-GENOME-SHA256-v1` raw/RC fingerprint audits
contradictory labels even when topology declarations differ. Exact or
reverse-complement-equivalent content:

- cannot occur with contradictory virus/host labels;
- fails by default when duplicated within one class; or
- can be deterministically reduced under `duplicate_policy = "drop"`, with the
  removed source and retained representative recorded.

Split planning groups equivalent content by that digest. This guards exact
content leakage, not all homology. Related but non-identical sequences require
the similarity protocol and appropriate external analysis.

Metadata is sequence-level and exhaustively joined: missing and unused rows are
errors. The ambiguous legacy date headings `deposited_at` and `create_date` are
rejected; only a verified first-public `release_date` is accepted. Topology is
declared per sequence. Grouped segments must agree in taxonomy. When every
segment has a release date, a grouped genome's public availability is the
latest of those dates. If any segment lacks a date, the whole group is treated
as unknown so partial metadata cannot make it appear historically available.
Additional metadata stays on each sequence; a group-level extra is retained
only when present with the same value on every segment. Users should retain
versioned accession identifiers because GenBank increments the version when a
sequence changes
([NCBI sequence identifiers](https://www.ncbi.nlm.nih.gov/genbank/sequenceids/)).

## Fragment sampling

For each source genome and requested length (L), eligible coordinates are all
0-based starts on contigs with at least (L) nucleotides. CHIMERA samples with
replacement from the union of those coordinates. A linear contig contributes
`contig_length - L + 1` coordinates; a declared circular contig contributes one
coordinate per source base because an interval may wrap the origin. Every
eligible start for a fixed length has equal probability. A fragment never
crosses between contigs.

The configured total `fragments_per_genome` is divided among requested lengths
as evenly as possible. When division has a remainder, a deterministic semantic
ranking selects which lengths receive one extra record. In `both` strand mode,
each accepted coordinate is deterministically assigned either its forward
substring or reverse complement; `forward` emits only the former.

Every truth interval is `[source_start, source_end)` on the forward source
contig. Linear ends are bounded by source length. Circular intervals use an
unwrapped end equal to start plus fragment length, so the end may exceed source
length when the interval wraps the declared origin. For a reverse-strand
record, the FASTA sequence is the reverse complement of the reconstructed
forward interval; the coordinates are not reversed. Fragment IDs are stable
opaque tokens and encode neither label nor source.

Candidates with a non-ACGT fraction above `max_ambiguous_fraction` are rejected.
Sampling has a finite attempt budget. Failure to obtain the requested records
raises an input error instead of silently shortening the dataset. Sampling
with replacement means identical coordinate or sequence observations are
possible; they remain separate generated observations with distinct IDs.

CHIMERA emits exact substrings. It does not reproduce abundance, community
structure, insert size, quality scores, sequencing errors, chimeric reads, or
paired-read correlation. CAMISIM provides community/read simulation with gold
standards ([DOI 10.1186/s13059-019-1593-1](https://doi.org/10.1186/s13059-019-1593-1)),
and InSilicoSeq models empirical short-read properties
([DOI 10.1093/bioinformatics/bty630](https://doi.org/10.1093/bioinformatics/bty630)).
Those are complementary references, not CHIMERA dependencies or validation
claims.

## Deterministic randomization

The master integer seed is never consumed as one global stateful stream.
CHIMERA derives domain-separated BLAKE2b sub-seeds from semantic parts such as
protocol, genome, length, fragment membership, and output order, then uses
Python `random.Random`. This prevents input iteration order and unrelated
protocol selection from shifting another decision. Ties are resolved with
stable IDs.

Input files and records are discovered and normalized in stable order; tables
and JSON use canonical ordering; gzip stores no source filename and has
`mtime=0`. These choices support repeatable semantic artifacts for identical
input bytes, resolved configuration, tool version, and compatible runtime.
They are not a promise that every future CHIMERA/Python implementation will
reproduce historical pseudorandom bytes. Archive the version, environment,
resolved configuration, input hashes, output checksums, and seed.

## Test 2A — random fragment split

1. Generate the complete requested fragment set from every genome.
2. For each genome independently, derive a stable shuffled order.
3. Choose a test count near `test_fraction`, bounded so train and test are both
   non-empty.
4. Semantically shuffle combined output so classes/sources are interleaved.

No fragment identifier occurs in both partitions. Every source genome occurs in
both partitions by design. This is a useful diagnostic for data loading,
optimization, and continuity with earlier evaluations, but it deliberately
permits genome-level dependence and must not be the sole generalization result.

## Test 2B — genome-level split

Genome content groups are partitioned before fragment generation. Viruses and
hosts are handled separately to approximate `test_fraction` while guaranteeing
at least one independent content group of each label in both partitions. Only
after the immutable plan is formed are fragments generated from the assigned
sources.

The validation invariant is zero overlap in fragment ID, user `genome_id`, and
canonical content hash. Test 2B answers whether a classifier recognizes unseen
source genomes in this sampled collection. It does not enforce a maximum
homology between different genomes.

## Test 2C — similarity-filtered split

### Candidate proposal

CHIMERA first makes the same kind of label-stratified, content-disjoint genome
proposal as 2B. Proposed training genomes remain training. Every proposed test
genome is a **candidate** and is compared with all training genomes. The
proposal is independent of input order and similarity values.

`candidate_test.fasta.gz` and its truth table contain the full proposal. This is
the denominator for similarity strata and exclusion accounting.

### Built-in MinHash screen

For each contig, CHIMERA enumerates unambiguous canonical k-mers; no k-mer spans
a contig boundary. Each canonical k-mer is hashed with stable BLAKE2b, and the
smallest `sketch_size` distinct hashes are retained. When both sketches retain
all unique k-mers, Jaccard is exact; otherwise consistent bottom-k sampling
estimates it.

For Jaccard (J > 0), CHIMERA reports the Mash-style estimate:

\[
D = -\frac{1}{k}\log\left(\frac{2J}{1+J}\right),
\qquad \widehat{I} = \operatorname{clip}(1-D, 0, 1).
\]

When (J=0), the result is “no detectable match,” not a measured identity of
zero. The lexically smallest training ID breaks an equal-similarity tie. The
method follows the MinHash/Mash rationale of Ondov et al.
([DOI 10.1186/s13059-016-0997-x](https://doi.org/10.1186/s13059-016-0997-x)),
but CHIMERA's implementation and defaults must be evaluated on the intended
data. Its estimate is not alignment ANI.

### External similarities and strict gate

For decisive publication analysis, an external table is an exact Cartesian
matrix with one row for every candidate-to-training pair. Every row carries
both genome IDs, `coverage_definition=aligned_fraction_shorter`, and a
versioned method identifier. A detected pair has fractional similarity and
aligned coverage. A nondetected pair retains both IDs and leaves both evidence
fields blank; a half-blank pair is invalid.

CHIMERA uses the greatest detected similarity for novelty assignment and
separately evaluates the strict compound gate against every detected pair. The
greatest qualifying pair supplies exclusion provenance even when another pair
is the numerical maximum. Lexical reference ID breaks ties. If every pair for
a candidate is nondetected, the result is `no_detectable_match` with no nearest
training ID. Thus no omitted pair can masquerade as a negative result, and all
candidates and training comparisons remain auditable.

A candidate is removed from strict test only when:

\[
similarity > max\_train\_similarity
\quad\text{and}\quad
(coverage\ is\ absent\ or\ coverage \ge min\_similarity\_coverage).
\]

Identity equality is retained; coverage equality qualifies. The built-in
screen has no alignment coverage, so its value is absent and the identity gate
applies. An external hit above identity but below minimum coverage remains in
strict test and carries an explicit assignment reason. Users must define why
that coverage logic is appropriate for their sequence type and external tool.

`test.*` contains the retained strict primary set. Above-gate exclusions remain
visible in `assignments.tsv` and `excluded.tsv`, and their fragments remain in
`candidate_test.*` and the relevant similarity stratum. Never report strict
accuracy without candidate/exclusion counts.

### Novelty strata

The selected training similarity described above assigns every candidate to one
mutually exclusive band. With defaults:

- `high_similarity`: `[0.90, 1]`;
- `moderate_similarity`: `[0.70, 0.90)`;
- `low_similarity`: `[0.30, 0.70)`;
- `distant_detectable`: measured `[0, 0.30)`;
- `no_detectable_match`: no detected similarity value.

These names are benchmark bins, not ranks or biological species definitions.
MIUViG discusses an operational vOTU boundary of at least 95% ANI over at least
85% of the shorter sequence
([DOI 10.1038/nbt.4306](https://doi.org/10.1038/nbt.4306)); that alignment-and-
coverage definition must not be conflated with a MinHash estimate, CHIMERA's
strict default, or a universal higher-rank threshold.

Recommended reporting includes performance and confidence intervals for every
non-empty candidate stratum, the strict set, excluded counts, external method
and version, reference snapshot, identity/coverage definitions, and the full
versioned pairwise table.

## Test 2D — temporal split

`release_date` is interpreted as first public accession release, consistent
with the NCBI Virus data report
([NCBI documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/virus/)).
For an explicit cutoff (t):

- train: effective release date `<= t`;
- test: effective release date `> t`;
- missing: error or explicit exclusion, never implicit training.

The cutoff is inclusive. If omitted, CHIMERA evaluates observed dates and
stably selects a class-viable cutoff nearest the requested test fraction. This
is a convenience for exploratory benchmarks; preregistered work should justify
an explicit date.

Exact content groups use the earliest represented public availability to keep a
later duplicate representation from leaking across the boundary. A segmented
genome's own date is the latest segment release because all grouped material is
needed to form that source.

### Retrospective caveat

Filtering today's reference collection by release date does not reconstruct a
historical database. Current accession sequence versions, taxonomic labels,
quality control, and which records survived curation can contain information
unavailable at (t). Therefore the standard CHIMERA output is accurately
described as a **release-date-filtered retrospective split**.

The prospective question—“would a model frozen at (t) recognize viruses
discovered later?”—requires the training reference, taxonomy, preprocessing,
and model to be frozen from an archived snapshot available at (t), with later
test material handled under a documented policy. CHIMERA can consume such
files but neither retrieves nor certifies their historical status. Report the
snapshot archive and identifiers before using prospective language.

## Test 2E — taxonomic holdout

CHIMERA chooses explicit or stably auto-selected viral values at one metadata
rank. All viral sources with those values test; remaining represented viral
taxa train. At least two viral values are required and at least one must remain
in training. Host sources are assigned independently using a label-stratified
genome/content-disjoint holdout so both labels are available in train and test.

Matching performs whitespace stripping and case-insensitive exact comparison
on the supplied rank value. It does not consult a taxonomy database, resolve
taxon IDs or ancestors, recognize synonyms, or merge spelling/name variants.
Users must normalize values against one named, versioned taxonomy snapshot and
audit polyphyly/reclassification risks. Broader environmental metadata can use
MIxS conventions ([canonical MIxS repository](https://github.com/GenomicsStandardsConsortium/mixs)).

2E measures generalization beyond represented string-defined taxon groups at
the chosen rank. It does not guarantee sequence dissimilarity; jointly inspect
2C similarity values. Conversely, a sequence-novel 2C item need not belong to a
new named taxon.

## Leakage and validation invariants

For every written protocol, CHIMERA checks fragment-ID disjointness and class
viability. For 2B–2E it additionally requires no train/test `genome_id` or
canonical source-content overlap. Protocol-specific plans validate temporal,
taxonomic, and similarity rules. The bundle stores these summaries in
`split.json` and an independent `chimera validate` command checks the completed
artifact and checksums. Validation reconstructs fragment sequence from the
per-sequence `sequences.tsv` and `source-sequences.fasta.gz` inventory, including
declared linear/circular topology and strand, and checks truth
`source_genome_id` and `source_content_group_id` against the retained sources.

The structured validation report separates primary train/test FASTA and truth
counts from auxiliary 2C `candidate_test`/stratum view counts. Its aggregate
record fields add those two categories and therefore count verified view rows,
not unique fragments or biological sources.

These controls address reference/assignment leakage inside the supplied
bundle. They do not detect:

- homologous or contaminated source sequence below the configured/external
  similarity analysis's sensitivity;
- duplicates in external pretraining corpora used by a model;
- taxonomy or date errors in source databases;
- hyperparameter selection on the final test set;
- features leaking through filenames, directory organization, side metadata,
  or downstream preprocessing;
- patient/site/time dependence not represented by `genome_id`.

Audit the model's entire training provenance, not only CHIMERA partitions.

## Reproducibility and reporting

The bundle records normalized references, source file hashes, resolved
configuration, a per-sequence source inventory and normalized source FASTA,
algorithm names/parameters, assignments, exclusions, truth, validation
summaries, tool version, seed, and file checksums. `execution.json` records
timestamps, Python, and platform and is intentionally excluded from semantic
checksums.

Public provenance is content-addressed: resolved configuration and receipts use
`sha256:<digest>` identifiers and semantic roles, not absolute input locations.
Creators should archive a controlled content-ID-to-accession/object mapping when
it is needed to retrieve original inputs.

The repository `environment.yml` gives bounded development dependencies, not a
fully resolved lock. Reproduction claims should archive an explicit lock/export
or immutable container digest in addition to the resolved CHIMERA configuration.

For a FAIR release, deposit immutable input references where licensing permits,
metadata, external similarity tables, the full bundle, environment lock or
container digest, and analysis code. Create a dataset-specific persistent
identifier and connect source/accession records with appropriate provenance.
The FAIR principles are described at
[DOI 10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18); DataCite's
current metadata schema documents version and related-identifier relations
([DataCite 4.7](https://schema.datacite.org/),
[relatedIdentifier](https://datacite-metadata-schema.readthedocs.io/en/4.7/properties/relatedidentifier/)).
CHIMERA does not currently emit DataCite, CAMI, or RO-Crate packages; converting
to those formats is an explicit downstream release step, not an implied
conformance claim.

## Validation status and limitations

CHIMERA's automated tests and bundle invariants verify implemented software
behavior. They do not constitute independent biological, clinical, regulatory,
or benchmark validation. Defaults are starting points, not universal choices.
Before publication, document why fragment lengths, ambiguity policy, split
fraction, similarity method/gates, coverage, temporal cutoff, taxonomy rank,
and source collection represent the target use case. Run sensitivity analyses
and preserve failures/exclusions.

Record those decisions in the deposited benchmark documentation and complete
the [dataset datasheet](../DATASHEET.md) with dataset-specific evidence.
