# Changelog

All notable user-visible changes are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and releases use
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

No user-visible changes recorded.

## [1.0.0] — 2026-08-08

This is the first publication-ready release of the rewritten CHIMERA engine.
It replaces the historical MetagenomeGenerator workflow with a schema-versioned,
leakage-aware benchmark bundle with descriptive protocol names and
purpose-first end-user documentation.

### Added

- One-command `chimera suite` orchestration of all five protocols, plus selected
  `chimera generate` runs and TOML/CLI override support.
- Strict recursive plain/gzip IUPAC DNA FASTA loading, mandatory metadata for
  multi-record FASTA, rejection of ambiguous legacy date columns, exhaustive
  sequence-level joins, multi-contig `genome_id` grouping, and normalized
  inventory inspection.
- Exact/reverse-complement-aware content hashing, cross-class conflict failure,
  and audited same-class duplicate handling.
- Deterministic uniform-coordinate fragment sampling with replacement,
  balanced requested lengths, both/forward strand modes, bounded ambiguity
  rejection, opaque label-free identifiers, and explicit linear/circular
  0-based half-open truth.
- Deterministic per-genome random-fragment diagnostic with non-empty
  partitions and interleaved semantic output order.
- Label-stratified whole-genome/content holdout.
- Genome-disjoint similarity candidate proposal, complete candidate view, strict
  identity/coverage gate, five candidate similarity strata, dependency-free
  canonical k-mer bottom-k MinHash/Mash-style estimates, and external
  all-candidate-vs-train similarity-table support.
- Inclusive first-public-release-date holdout with explicit missing-data
  accounting, stable viable auto-cutoff, and retrospective provenance wording.
- Explicit or stable automatic viral taxon holdout with independently
  genome-disjoint host partitions.
- Atomic bundle publication and guarded `--force`, deterministic JSON/TSV/gzip,
  semantic BLAKE2b seed derivation, content-addressed input receipts without
  absolute local paths, resolved configuration, per-sequence source inventory,
  assignments, exclusions, manifests, truth, embedded schemas, checksums, and
  human report. Manifests include a software-content SHA-256 receipt that stays
  available for installed wheels even when Git metadata is absent.
- Independent `chimera validate`, `chimera inspect`, and table/JSON
  `chimera schema` discovery commands with documented exit statuses.
- Typed Python package, automated tests and coverage threshold, Ruff, strict
  mypy, required CI across Python 3.11–3.14, wheel/sdist inspection,
  pre-commit, Make targets, bounded Conda development environment, non-root OCI
  container, and
  security policy.
- Scientific methodology and user documentation, dataset datasheet,
  citation/CodeMeta records, contributor governance, and a runnable network-free
  tiny fixture.

### Changed

- The maintained executable is now `chimera`.
- Generated FASTA headers no longer expose source or class; downstream tools
  must join opaque IDs to truth tables.
- Split assignment is explicit at the biological unit required by each claim,
  rather than relying on one implicit random record split.
- Output is a versioned directory bundle rather than a legacy flat/ad-hoc
  layout. Existing unrecognized directories are never replaced.

### Deprecated

- `metagenome-generator` remains as a compatibility console entry point and
  prints a deprecation warning before delegating to the new CLI. It may be
  removed in a future major release.

### Known limitations

- Generated records are exact synthetic substrings, not empirical read/error,
  paired-end, abundance, or community simulations.
- Built-in similarity is a Mash-style MinHash identity estimate, not
  alignment-derived ANI; a versioned external table is recommended for
  publication analyses requiring ANI/coverage semantics.
- Temporal filtering of current references is retrospective unless the user
  supplies an authentic immutable historical reference and taxonomy snapshot.
- Taxonomic matching is case-insensitive exact-string comparison and does not
  resolve taxonomy IDs, lineages, synonyms, or renamed taxa.
- CHIMERA validates bundle mechanics and declared invariants, not the biological
  correctness of user-supplied labels, dates, taxonomy, or model pretraining
  provenance.

[Unreleased]: https://github.com/Alexander-Mitrofanov/CHIMERA/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/Alexander-Mitrofanov/CHIMERA/tree/v1.0.0
