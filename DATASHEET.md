# Datasheet for CHIMERA-generated benchmark datasets

This document is the project-level datasheet template for bundle schema
`urn:chimera:benchmark-bundle:1`. It describes what CHIMERA generates and the
bundled tiny fixture. A research release must copy this file into its deposited
dataset and replace every bracketed field with dataset-specific facts.

The structure is informed by *Datasheets for Datasets*
([DOI 10.1145/3458723](https://doi.org/10.1145/3458723)); inclusion here does
not imply review or endorsement by that work's authors.

## Dataset identity

- **Name:** [dataset title]
- **Version:** [immutable dataset version]
- **Persistent identifier:** [dataset DOI/handle; do not insert the CHIMERA
  software URL or a methods-paper DOI]
- **Bundle schema:** `urn:chimera:benchmark-bundle:1`
- **CHIMERA version:** [for example `1.0.0`]
- **Creators and affiliations:** [names/ORCIDs/roles]
- **Maintainer:** [durable contact or repository]
- **Release date:** [YYYY-MM-DD]
- **License/access terms:** [data-specific license and any source restrictions]
- **Source snapshots:** [database releases, retrieval dates, accession lists,
  URLs/API commands, and persistent identifiers]
- **Software citation:** See `CITATION.cff`; cite software separately.

CHIMERA 1.0.0 has no asserted project DOI. Each published bundle should receive
its own persistent identifier where possible.

## Motivation

CHIMERA bundles are intended to evaluate binary virus-versus-host fragment
classifiers under five explicitly different train/test relationships:

1. random fragments from seen genomes as a diagnostic;
2. entirely held-out source genomes;
3. genome-disjoint candidates measured and filtered by training similarity;
4. first-public-release dates on opposite sides of a cutoff; and
5. selected viral taxonomy values absent from training.

The suite exists because a random fragment split cannot alone support an
unseen-genome or novel-virus claim. The intended deliverable is a transparent
benchmark with truth, assignments, exclusions, resolved parameters, and
checksums—not a trained model or a universal biological reference standard.

Dataset-specific motivation: [state the concrete biological question,
deployment context, why virus/host classes and these reference populations are
appropriate, and the primary protocol/endpoint].

## Composition

A complete bundle contains retained genome and per-sequence inventories, the
normalized `source-sequences.fasta.gz`, preflight exclusions, five protocol
directories, compressed fragment FASTA, one-to-one fragment truth, source-level
assignment/exclusion tables, split and bundle manifests, a human report,
resolved configuration, an embedded schema snapshot, execution facts, and
SHA-256 checksums. See
[`docs/OUTPUT_FORMATS.md`](docs/OUTPUT_FORMATS.md).

Each FASTA record is an exact synthetic DNA substring of one supplied source
contig, optionally reverse-complemented. It has:

- one binary ground-truth label, `virus` or `host`;
- an opaque label-free fragment ID;
- exact `source_genome_id`, canonical `source_content_group_id`, source sequence
  ID, orientation, and linear/circular 0-based half-open coordinates in truth;
- an exact requested nucleotide length; and
- a semantic `partition`, serialized `view`, and protocol-specific similarity
  metadata.

These records are not empirical reads. They have no qualities, abundance,
paired-end relationship, empirical error process, or simulated community
profile.

Dataset-specific composition:

- Retained sources: [N viruses; N hosts; taxonomic and length distributions]
- Preflight exclusions: [count/reasons]
- Fragment lengths/count per source: [values]
- Random-fragment train/test records and sources: [counts]
- Genome-holdout train/test records and sources: [counts]
- Similarity candidate/strict/excluded/stratum counts: [counts]
- Temporal train/test/excluded counts and cutoff: [counts/date]
- Taxonomic train/test/excluded counts, rank, and values: [counts/values]
- Missing values: [fields/counts/policy]
- Sensitive attributes: [assessment]

Use `manifest.json`, `references.tsv`, each `split.json`, and exclusion tables
to populate these counts. Do not count only FASTA records; fragments sharing a
source genome are biologically dependent.

## Data sources and acquisition

CHIMERA performs no download. Bundle creators supply all FASTA and metadata and
are responsible for accuracy, consent/access, licensing, and snapshot
provenance. Document for every upstream collection:

- database/provider and release or immutable snapshot;
- retrieval date, exact query/API/command, and raw file hash;
- original license/terms and redistribution decision;
- inclusion/exclusion/quality criteria;
- accession.version values and mapping to CHIMERA IDs;
- how virus and host/non-viral labels were established;
- taxonomy database/version and normalization procedure;
- how `release_date` was obtained and distinguished from collection date; and
- any contamination, completeness, assembly, or host-association screening.

NCBI states that Virus `releaseDate` is first public release
([data report documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/virus/)).
GenBank accession versions change when sequence content changes
([identifier documentation](https://www.ncbi.nlm.nih.gov/genbank/sequenceids/)).
Pin both content and identifiers.

Dataset-specific acquisition record: [complete the items above].

## Preprocessing and labeling

CHIMERA recursively discovers supported plain/gzip DNA FASTA, validates unique
stable IDs, removes whitespace, uppercases the ungapped IUPAC DNA alphabet, and
rejects RNA `U`, gaps, invalid symbols, empty records, and malformed headers.
A multi-record FASTA requires explicit metadata; CHIMERA does not infer whether
its records are independent or segmented. Metadata is exhaustively joined by
`sequence_id`; contigs sharing `genome_id` form one source. `topology` is
sequence-level, and taxonomy must be consistent within a source. Ambiguous
legacy date columns `deposited_at` and `create_date` are rejected in favor of a
verified first-public `release_date`.

Topology-aware, reverse-complement-invariant digest v2 defines source grouping;
circular origins are rotation invariant, while linear/circular declarations are
distinct. A separate topology-agnostic exact raw/RC fingerprint makes
cross-class equivalence fatal even when topology metadata differs. Same-class
duplicates under the same topology semantics fail or are deterministically
excluded and reported. These checks do not establish that a reference is
biologically uncontaminated or correctly labeled.

The input channel (`--virus` or `--host`) defines the class; an optional
metadata label is a consistency assertion. Label quality is therefore inherited
from the creator's source selection and upstream databases.

Dataset-specific preprocessing beyond CHIMERA: [adapter/contamination
filtering, sequence selection, deduplication before CHIMERA, taxonomy mapping,
identifier pseudonymization, scripts/versions, and manual decisions].

## Generation and split construction

Fragments are sampled with replacement from uniformly weighted eligible start
coordinates at each requested length; lengths are approximately balanced per
genome, fragments never cross contigs, and ambiguous candidates undergo bounded
rejection. Linear intervals end within the source; declared circular intervals
may wrap the origin and use an unwrapped end coordinate. Both orientations are
sampled by default. Domain-separated BLAKE2b sub-seeds and stable ordering
isolate decisions from input order.

### Random-fragment diagnostic

Every genome contributes non-empty train and test fragments. Fragment IDs are
disjoint; source overlap is deliberate. This protocol must not be presented as
evidence for unseen-genome generalization.

### Genome holdout

Whole user/content groups are assigned label-stratified before fragment
generation. Source ID and exact canonical content do not cross partitions.

### Similarity-filtered holdout

A content-disjoint candidate proposal is compared with training. The complete
proposal is retained as `candidate_test`; `test` is the primary strict set.
Candidate similarity strata are `high_similarity`, `moderate_similarity`,
`low_similarity`, `distant_detectable`, and `no_detectable_match`. Strict and
candidate sets must be reported separately.

The built-in canonical k-mer MinHash reports a Mash-style identity estimate,
not alignment ANI. Dataset creators should archive a versioned external
alignment-derived all-candidate-vs-train table for publication analyses where
ANI/coverage is the intended construct. Record tool/version, inputs, commands,
coverage definition, no-match policy, and thresholds.

### Temporal holdout

Training dates are on/before an inclusive cutoff; test dates are later. In the
absence of an actual immutable historical reference and taxonomy snapshot, this
is a **release-date-filtered retrospective** split, not prospective discovery.

### Taxonomic holdout

Selected supplied viral rank strings are excluded from training. Matching is
case-insensitive but otherwise exact and does not resolve taxonomy IDs,
lineages, synonyms, or renamed values. The dataset creator must normalize a
versioned taxonomy. Taxonomic novelty is not evidence of sequence novelty.

Dataset-specific generation record: [attach exact TOML, CLI invocation,
standard output/error, seed, CHIMERA/Python versions, dependency lock/container
digest, and sensitivity analyses].

The repository `environment.yml` is a bounded contributor environment, not a
fully solved lockfile. Archive the actual resolved lock/export or immutable
container digest used for the dataset.

## Quality assurance

CHIMERA preflights source/class counts and length feasibility, validates each
plan, checks primary partition fragment IDs, source IDs, exact content leakage,
and class viability, writes atomic deterministic files, and hashes semantic
bundle contents. `chimera validate` independently checks a completed bundle,
reconstructing fragments against the normalized per-sequence inventory and
source FASTA. Its structured report distinguishes primary train/test
FASTA/truth counts from auxiliary similarity candidate/stratum view counts; aggregate
record counts include both categories and are not unique-source counts.

Dataset release evidence:

- [ ] `chimera validate BUNDLE --json` report attached
- [ ] `checksums.sha256` verified after deposition
- [ ] independent clean-environment rerun compared
- [ ] FASTA/truth coordinate reconstruction sampled or exhaustive
- [ ] external similarity table coverage/completeness checked
- [ ] excluded/missing records manually reviewed
- [ ] source label/taxonomy/date spot audit recorded

Software tests verify implementation properties, not biological truth or
fitness for a clinical/regulatory purpose.

## Intended uses

Appropriate uses, subject to source suitability, include:

- controlled comparison of virus/host fragment classifiers;
- measuring degradation from seen fragments toward unseen, dissimilar, later,
  or held-out-taxonomy sources;
- data-loader and leakage diagnostic work with random-fragment splits;
- ablation/sensitivity analysis across fragment lengths and similarity strata; and
- auditable reproduction of a specifically configured published evaluation.

Dataset-specific intended users and decisions: [describe].

## Out-of-scope and discouraged uses

Do not use a CHIMERA bundle by itself to:

- claim clinical, diagnostic, biosafety, or regulatory validity;
- estimate real-world prevalence or community composition;
- represent empirical platform errors/read qualities or complete metagenomes;
- claim prospective temporal discovery without a genuine historical snapshot;
- infer taxonomic rank from MinHash bins or treat defaults as universal viral
  boundaries;
- report random-fragment performance as unseen-genome performance;
- infer organism pathogenicity, host association, phenotype, or ecological
  risk from the binary label;
- redistribute source sequences contrary to their terms; or
- benchmark a model whose pretraining/reference corpus overlaps test sources
  without disclosure and appropriate controls.

Users must evaluate whether the source population, class definition, fragment
process, and partitions represent their deployment setting.

## Known biases, risks, and limitations

- Reference databases overrepresent cultivable, clinically studied, and
  geographically/institutionally sampled organisms and underrepresent unknown
  diversity.
- Database virus/host labels, dates, assemblies, and taxonomy can be wrong or
  revised. Exact accession/version and snapshot provenance are essential.
- Equal fragments per genome do not model environmental abundance and can
  overweight short sources or repeated coordinates.
- Fragments from one genome are dependent; fragment-level confidence intervals
  may be anticonservative.
- Exact-content grouping does not remove all homology. Built-in MinHash has
  k-mer/sketch/fragmentation/ambiguity sensitivity and is not alignment ANI.
- External “best hit” behavior inherits its upstream method's definition,
  coverage, database, and failure modes.
- Release-date filtering of a current database contains survivorship and
  present-day annotation information.
- Exact supplied taxonomy strings can split synonyms or merge improperly
  normalized names; rank labels can be unstable.
- Binary `host` collapses diverse non-viral content and may not represent the
  target negative distribution.
- Opaque generated IDs reduce direct label leakage but cannot prevent leakage
  through downstream paths, metadata, pretraining, or feature engineering.
- Small or empty similarity strata limit inference. Exclusions change the strict
  target population and must be disclosed.

Dataset-specific bias/risk assessment and mitigations: [describe].

## Ethical, legal, and security considerations

Public genome references are not automatically free of legal, ethical, dual-use,
or privacy concerns. Environmental/clinical metadata can reveal locations,
patients, facilities, outbreaks, or protected traditional knowledge even when
sequences appear non-identifying. Assess consent, access agreements,
benefit-sharing, indigenous data sovereignty, export controls, embargoes, and
dual-use review for the actual sources and intended release.

CHIMERA does not write absolute local input paths into the bundle. It replaces
them with `sha256:<digest>` content IDs and semantic roles in the resolved
configuration, manifest receipts, and source inventories. Preserve a controlled
immutable mapping from those IDs to accessions or source objects when retrieval
is required. Content hashes and user-supplied extra metadata can still enable
correlation or reveal sensitive context, so review the complete bundle before
deposit.

Dataset-specific assessment, review body, and restrictions: [describe].

## Distribution, licensing, and access

CHIMERA software is MIT-licensed; that license does **not** relicense input
references or generated sequence-derived data. The dataset creator must select
accurate data-specific terms and propagate upstream attribution/restrictions.

- Repository/deposit: [URL/identifier]
- Access procedure: [open/registered/controlled]
- Dataset license: [SPDX or full terms]
- Upstream restrictions: [list]
- Takedown/correction process: [contact/process]
- Retention policy: [describe]

Where conversion to CAMI, DataCite, MIxS, or RO-Crate is performed downstream,
name the exact version and validator. The native CHIMERA bundle does not claim
conformance to those packaging standards.

## Maintenance and versioning

The dataset maintainer—not the CHIMERA software project by default—is
responsible for corrections, source-database updates, and user questions.

- Maintenance period: [dates/commitment]
- Contact: [durable address/repository]
- Issue and correction policy: [describe]
- Version semantics: [describe]
- Deprecation/takedown policy: [describe]

Never mutate a deposited bundle in place. A reference, metadata, similarity
table, configuration, algorithm, or correction change produces a new dataset
version, new checksums, and updated datasheet/changelog. Link versions and
sources using accurate persistent-identifier relations.

## Tiny fixture declaration

`examples/tiny/` is a software smoke-test/demo, not a scientific benchmark. It
contains four deterministic synthetic 320 nt strings: two labeled viruses in
two fictional families and two labeled hosts. Each label has one release date
on either side of `2020-12-31`; all accessions, taxa, species, dates, and
sequences are invented. It emits eight records per genome at 31 nt and 61 nt.
It may be redistributed under the repository license and must not be cited as
biological evidence.
