# Contributing to CHIMERA

Thank you for helping make CHIMERA scientifically reliable and easier to use.
Contributions include code, tests, documentation, benchmark fixtures, metadata
schema review, reproducibility reports, and careful bug reports.

Report security-sensitive issues through [SECURITY.md](SECURITY.md), not a
public issue.

## Before starting

For a defect, search existing issues and provide a minimal synthetic reproducer.
For a new flag, split protocol, schema change, dependency, or scientifically
meaningful default, open a design issue first. Describe:

- the biological/evaluation claim the change supports;
- the assignment unit and leakage threat model;
- proposed CLI/config/schema compatibility;
- determinism and provenance implications;
- failure behavior and missing-data policy;
- evidence or primary references supporting the method; and
- migration and release impact.

Do not attach controlled, embargoed, identifiable, or license-restricted genome
data. Minimize real biological sequence in reports; a synthetic fixture is
preferred.

## Development setup

CHIMERA supports the Python versions declared in `pyproject.toml`; its small
runtime dependency set is declared there. Create an isolated environment with
Python 3.11 or newer:

```console
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
make install-dev
python -m pre_commit install
```

The conda environment is an alternative for contributors:

```console
conda env create -f environment.yml
conda activate chimera-dev
```

## Make targets

```console
make help          # list targets
make format        # apply Ruff lint fixes and formatting
make lint          # configured Ruff rules
make typecheck     # strict mypy
make test          # complete pytest suite
make coverage      # branch coverage and configured threshold
make ci            # stable local CI gates plus package inspection
make smoke         # build, inspect, install wheel, and run CLI smoke tests
make security      # audit installed dependencies
make container     # build the non-root OCI image
make pre-commit    # all hooks against all files
```

Run the narrowest relevant tests during iteration, then `make ci` before a pull
request. A release-affecting change should also pass `make smoke` and the tiny
offline workflow:

```console
chimera suite --config examples/tiny/chimera.toml --outdir /tmp/chimera-tiny-check
chimera validate /tmp/chimera-tiny-check
```

Choose a new temporary output name or a safe temporary directory rather than
using `--force` against valuable data.

## Engineering expectations

- Support Python 3.11+ and keep runtime dependencies minimal, justified, and
  narrowly bounded.
- Use strict type annotations and concise public docstrings.
- Prefer immutable models, explicit errors, bounded work, stable semantic IDs,
  and atomic writes.
- Derive random decisions from domain-separated semantic seeds. Do not depend
  on input order, process-randomized `hash()`, wall-clock time, or global PRNG
  state.
- Preserve user data. Never overwrite an unrecognized output directory or
  silently omit requested records, metadata failures, or exclusions.
- Keep generated FASTA headers opaque. Labels and source provenance belong in
  truth tables, never an identifier or directory convention visible to a model.
- Treat paths and FASTA/TSV/TOML/JSON as untrusted input; include file/line
  context in actionable errors where possible.
- Use standard-library implementations where clear and maintainable. If an
  external scientific algorithm is needed, make method/version/provenance
  explicit and avoid shell interpolation.
- Keep user-facing flags intuitive, stable, and consistent between TOML and CLI.
  New defaults are scientific API changes and require justification.

## Scientific-method changes

A split, similarity, coordinate, metadata, or sampling change requires more
than unit coverage. Include:

1. a written invariant and threat model;
2. focused tests, boundary cases, order/permutation tests, and negative tests;
3. a deterministic golden or property-based check where appropriate;
4. validation of class viability, overlap, exclusions, and failure behavior;
5. updates to methodology, user guide, output schema, dataset datasheet, and
   changelog as applicable;
6. primary-source citations without implying external endorsement; and
7. a migration/schema-version decision for existing bundles and consumers.

Defaults must not be presented as universal biological thresholds. In
particular, MinHash similarity is not alignment ANI, taxonomy strings require a
versioned normalization source, and present-day release-date filtering is not a
prospective historical experiment.

## Tests and fixtures

Tests should be fast, deterministic, network-free, and isolated with temporary
directories. Use clearly fictional IDs, taxonomy, dates, and synthetic IUPAC
DNA. Do not use a real pathogen sequence merely for convenience.

At minimum, cover:

- happy path and exact error message/context;
- minimum/maximum and equality boundaries;
- input record/file permutations;
- both labels and multi-contig sources;
- ambiguity, reverse complement, duplicate/conflict, and missing metadata;
- train/test ID, source ID, and content-hash invariants;
- rerun byte stability of semantic files; and
- atomic failure/no-partial-output behavior.

Do not weaken the configured branch-coverage threshold to land a change.

## Documentation and metadata

Commands in documentation must run from a clean checkout or clearly state
their working directory. Name exact schemas, operators, coordinate systems,
defaults, caveats, and output views. Check external links and prefer DOI links
or primary official documentation.

When releasing, synchronize version/date and scope across `pyproject.toml`,
`src/chimera/__init__.py`, `CHANGELOG.md`, `CITATION.cff`, and `codemeta.json`.
Never invent a DOI, ORCID, affiliation, funding source, validation status, or
author identity. Dataset citations are separate from the software citation.

## Pull requests

Keep changes focused and avoid unrelated formatting. A pull request should
include:

- a concise problem/outcome statement;
- linked issue or design decision where relevant;
- tests and commands run;
- user/scientific/reproducibility impact;
- schema/CLI/migration impact;
- documentation and changelog updates; and
- confirmation that fixtures contain no restricted or sensitive material.

Reviewers may ask for additional domain review when a change affects biological
interpretation. Approval means the change meets project standards; it does not
certify clinical or universal biological validity.

## Reporting bugs

Include CHIMERA/Python/OS versions, installation method, exact command and
configuration (with private paths/data removed), full diagnostic text, expected
and observed behavior, and a minimal synthetic input. For nondeterminism, attach
both `resolved-config.json` files, manifests, checksums, and which semantic files
differ. For scientific concerns, state the expected invariant and supporting
reference.

## License

By submitting a contribution, you agree that it is licensed under the
repository's [MIT License](LICENSE). You must have the right to contribute all
included code, documentation, and data. Retain required attribution for adapted
material and identify its license.
