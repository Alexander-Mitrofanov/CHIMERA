"""Small, stable command-line interface for CHIMERA."""

from __future__ import annotations

import argparse
import json
import logging
import sys
import tomllib
from collections.abc import Sequence
from pathlib import Path
from typing import Any, NoReturn

from . import __version__
from .config import BenchmarkConfig
from .errors import ChimeraError, ConfigurationError, IntegrityError

LOGGER = logging.getLogger("chimera")


class _ArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> NoReturn:
        self.print_usage(sys.stderr)
        self.exit(2, f"chimera: error: {message}\n")


def _path(value: str) -> Path:
    return Path(value).expanduser()


def _fraction(value: str) -> float:
    try:
        result = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number between 0 and 1") from exc
    if not 0.0 <= result <= 1.0:
        raise argparse.ArgumentTypeError("must be a number between 0 and 1")
    return result


def _positive_int(value: str) -> int:
    try:
        result = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if result < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return result


def _add_generation_arguments(
    parser: argparse.ArgumentParser, *, allow_split_selection: bool
) -> None:
    input_group = parser.add_argument_group("inputs")
    input_group.add_argument(
        "--virus",
        dest="virus_paths",
        action="append",
        type=_path,
        metavar="FASTA_OR_DIR",
        default=argparse.SUPPRESS,
        help="Viral FASTA file/directory; repeat for multiple inputs.",
    )
    input_group.add_argument(
        "--host",
        dest="host_paths",
        action="append",
        type=_path,
        metavar="FASTA_OR_DIR",
        default=argparse.SUPPRESS,
        help="Host/non-viral FASTA file/directory; repeat for multiple inputs.",
    )
    input_group.add_argument(
        "--metadata",
        dest="metadata_path",
        type=_path,
        metavar="TSV",
        default=argparse.SUPPRESS,
        help="Sequence metadata TSV (required for temporal/taxonomic protocols).",
    )
    input_group.add_argument(
        "--config",
        type=_path,
        metavar="TOML",
        help="TOML configuration; explicit flags override its values.",
    )

    output_group = parser.add_argument_group("output and reproducibility")
    output_group.add_argument(
        "--outdir",
        dest="output_dir",
        type=_path,
        metavar="DIR",
        default=argparse.SUPPRESS,
        help="New benchmark-bundle directory (default: benchmark-output).",
    )
    output_group.add_argument(
        "--seed",
        type=int,
        default=argparse.SUPPRESS,
        help="Master deterministic seed (default: 42).",
    )
    output_group.add_argument(
        "--force",
        dest="overwrite",
        action="store_true",
        default=argparse.SUPPRESS,
        help="Replace an existing outdir using an atomic backup/commit.",
    )
    output_group.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs and resolve all split plans without writing fragments.",
    )

    protocol_group = parser.add_argument_group("evaluation protocols")
    if allow_split_selection:
        protocol_group.add_argument(
            "--split",
            dest="splits",
            action="append",
            metavar="NAME",
            default=argparse.SUPPRESS,
            help=(
                "Protocol: random, genome, similarity, temporal, taxonomy, or all; "
                "repeat or comma-separate (default: all)."
            ),
        )
    protocol_group.add_argument(
        "--test-fraction",
        type=_fraction,
        default=argparse.SUPPRESS,
        help="Target genome/fragment test fraction (default: 0.20).",
    )
    protocol_group.add_argument(
        "--release-date-cutoff",
        dest="temporal_cutoff",
        metavar="YYYY-MM-DD",
        default=argparse.SUPPRESS,
        help="Inclusive training cutoff; by default CHIMERA selects a viable cutoff.",
    )
    protocol_group.add_argument(
        "--holdout-rank",
        dest="taxonomy_rank",
        metavar="RANK",
        default=argparse.SUPPRESS,
        help="Taxonomic rank to exclude from training (default: family).",
    )
    protocol_group.add_argument(
        "--holdout-taxon",
        dest="holdout_taxa",
        action="append",
        metavar="NAME",
        default=argparse.SUPPRESS,
        help="Taxon to hold out; repeat. Stable auto-selection is used if omitted.",
    )
    protocol_group.add_argument(
        "--auto-holdout-count",
        type=_positive_int,
        default=argparse.SUPPRESS,
        help="Number of viral taxa to auto-hold out (default: 1).",
    )
    protocol_group.add_argument(
        "--missing-metadata",
        choices=("error", "exclude"),
        default=argparse.SUPPRESS,
        help="Fail or explicitly exclude records lacking required split metadata.",
    )

    fragment_group = parser.add_argument_group("fragment generation")
    fragment_group.add_argument(
        "--fragment-length",
        dest="fragment_lengths",
        action="append",
        type=_positive_int,
        metavar="NT",
        default=argparse.SUPPRESS,
        help="Exact fragment length in nucleotides; repeat (default: 500).",
    )
    fragment_group.add_argument(
        "--fragments-per-genome",
        type=_positive_int,
        metavar="N",
        default=argparse.SUPPRESS,
        help="Total fragments emitted per source genome (default: 100).",
    )
    fragment_group.add_argument(
        "--strand",
        dest="strand_mode",
        choices=("both", "forward"),
        default=argparse.SUPPRESS,
        help="Sample both orientations or forward only (default: both).",
    )
    fragment_group.add_argument(
        "--max-ambiguous-fraction",
        type=_fraction,
        metavar="FRACTION",
        default=argparse.SUPPRESS,
        help="Maximum non-ACGT fraction in an emitted fragment (default: 0.05).",
    )
    fragment_group.add_argument(
        "--duplicate-policy",
        choices=("error", "drop"),
        default=argparse.SUPPRESS,
        help="Same-class genome-content duplicate policy (default: error).",
    )

    similarity_group = parser.add_argument_group("similarity protocol")
    similarity_group.add_argument(
        "--similarity-table",
        type=_path,
        metavar="TSV",
        default=argparse.SUPPRESS,
        help="External all-test-vs-train similarities on [0,1] (e.g. skani/FastANI).",
    )
    similarity_group.add_argument(
        "--similarity-k",
        type=_positive_int,
        metavar="K",
        default=argparse.SUPPRESS,
        help="Odd canonical k-mer length for built-in MinHash (default: 21).",
    )
    similarity_group.add_argument(
        "--sketch-size",
        type=_positive_int,
        metavar="N",
        default=argparse.SUPPRESS,
        help="Bottom-k sketch size (default: 2000).",
    )
    similarity_group.add_argument(
        "--max-train-similarity",
        type=_fraction,
        metavar="FRACTION",
        default=argparse.SUPPRESS,
        help="Strict test maximum; higher candidates are accounted as excluded (default: .95).",
    )
    similarity_group.add_argument(
        "--min-similarity-coverage",
        type=_fraction,
        metavar="FRACTION",
        default=argparse.SUPPRESS,
        help=(
            "Minimum aligned fraction for an external-table hit to trigger exclusion "
            "(default: .85)."
        ),
    )
    similarity_group.add_argument(
        "--similarity-high",
        type=_fraction,
        default=argparse.SUPPRESS,
        help="Lower boundary of the high-similarity stratum (default: .90).",
    )
    similarity_group.add_argument(
        "--similarity-moderate",
        type=_fraction,
        default=argparse.SUPPRESS,
        help="Lower boundary of the moderate stratum (default: .70).",
    )
    similarity_group.add_argument(
        "--similarity-low",
        type=_fraction,
        default=argparse.SUPPRESS,
        help="Detection boundary of the low stratum (default: .30).",
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the public parser without performing I/O."""

    parser = _ArgumentParser(
        prog="chimera",
        description=(
            "Generate leakage-aware virus-versus-host benchmark datasets with auditable truth."
        ),
    )
    parser.add_argument("--version", action="version", version=f"CHIMERA {__version__}")
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
        help="Diagnostic verbosity on stderr (default: INFO).",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress progress messages.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    suite = subparsers.add_parser(
        "suite",
        help="Generate all five benchmark protocols in one reproducible bundle (recommended).",
        description="Generate the full leakage-aware benchmark suite in one command.",
    )
    _add_generation_arguments(suite, allow_split_selection=False)

    generate = subparsers.add_parser(
        "generate",
        help="Generate one or more selected evaluation protocols.",
        description="Generate selected leakage-aware evaluation datasets.",
    )
    _add_generation_arguments(generate, allow_split_selection=True)

    validate = subparsers.add_parser(
        "validate", help="Independently validate a generated benchmark bundle."
    )
    validate.add_argument("bundle", type=_path, metavar="DIR")
    validate.add_argument(
        "--json", action="store_true", help="Print the validation report as JSON."
    )

    inspect = subparsers.add_parser(
        "inspect", help="Preflight references and print their normalized inventory."
    )
    inspect.add_argument("--virus", dest="virus_paths", action="append", required=True, type=_path)
    inspect.add_argument("--host", dest="host_paths", action="append", required=True, type=_path)
    inspect.add_argument("--metadata", dest="metadata_path", type=_path)
    inspect.add_argument("--duplicate-policy", choices=("error", "drop"), default="error")
    inspect.add_argument("--json", action="store_true")

    from .schema_resources import JSON_SCHEMA_NAMES

    schema = subparsers.add_parser(
        "schema", help="Print a TSV header contract or versioned JSON Schema."
    )
    schema.add_argument(
        "name",
        choices=("metadata", "truth", "references", *JSON_SCHEMA_NAMES),
        nargs="?",
        default="metadata",
    )
    return parser


def _load_toml_mapping(path: Path) -> dict[str, Any]:
    try:
        with path.open("rb") as handle:
            document = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise ConfigurationError(f"Cannot read TOML configuration {path}: {exc}") from exc
    section = document.get("benchmark", document)
    if not isinstance(section, dict):
        raise ConfigurationError("TOML [benchmark] must be a table")
    return dict(section)


def _config_from_namespace(args: argparse.Namespace) -> BenchmarkConfig:
    values: dict[str, Any] = {}
    config_path: Path | None = getattr(args, "config", None)
    base_dir = Path.cwd()
    if config_path is not None:
        config_path = config_path.resolve()
        values.update(_load_toml_mapping(config_path))
        base_dir = config_path.parent
    cli_values = vars(args).copy()
    for key in (
        "command",
        "config",
        "dry_run",
        "log_level",
        "quiet",
        "similarity_high",
        "similarity_moderate",
        "similarity_low",
    ):
        cli_values.pop(key, None)
    # CLI paths are relative to the invocation directory, not to a config file.
    for key in ("virus_paths", "host_paths"):
        if key in cli_values:
            cli_values[key] = [path.resolve() for path in cli_values[key]]
    for key in ("metadata_path", "similarity_table", "output_dir"):
        if key in cli_values:
            cli_values[key] = cli_values[key].resolve()
    values.update(cli_values)
    band_keys = {
        "high": getattr(args, "similarity_high", None),
        "moderate": getattr(args, "similarity_moderate", None),
        "low": getattr(args, "similarity_low", None),
    }
    if any(value is not None for value in band_keys.values()):
        existing = values.get("similarity_bands", {})
        if not isinstance(existing, dict):
            raise ConfigurationError("similarity_bands in TOML must be a table")
        values["similarity_bands"] = {
            **existing,
            **{key: value for key, value in band_keys.items() if value is not None},
        }
    if args.command == "suite":
        requested_splits = values.get("splits", "all")
        from .config import SplitKind, parse_splits

        if parse_splits(requested_splits) != tuple(SplitKind):
            raise ConfigurationError(
                "chimera suite always generates all five protocols; use 'chimera generate' "
                "for a selected subset"
            )
        values["splits"] = "all"
    return BenchmarkConfig.from_mapping(values, base_dir=base_dir)


def _run_generate(args: argparse.Namespace) -> int:
    from .pipeline import generate_benchmark

    config = _config_from_namespace(args)
    result = generate_benchmark(config, dry_run=args.dry_run)
    if not args.quiet:
        if args.dry_run:
            LOGGER.info("Preflight passed: %s", result.summary)
        else:
            LOGGER.info("Benchmark bundle ready: %s", result.output_dir)
    print(json.dumps(result.as_dict(), sort_keys=True))
    return 0


def _run_validate(args: argparse.Namespace) -> int:
    from .validation import validate_bundle

    report = validate_bundle(args.bundle)
    if args.json:
        print(json.dumps(report.as_dict(), indent=2, sort_keys=True))
    else:
        print(f"PASS: {report.summary()}")
    return 0


def _run_inspect(args: argparse.Namespace) -> int:
    from .config import DuplicatePolicy
    from .output import reference_rows
    from .references import load_reference_catalog

    catalog = load_reference_catalog(
        tuple(path.resolve() for path in args.virus_paths),
        tuple(path.resolve() for path in args.host_paths),
        metadata_path=args.metadata_path.resolve() if args.metadata_path else None,
        duplicate_policy=DuplicatePolicy(args.duplicate_policy),
    )
    rows = reference_rows(catalog.genomes)
    if args.json:
        print(json.dumps({"references": rows}, indent=2, sort_keys=True))
    else:
        print("genome_id\tlabel\tlength_nt\tcontigs\trelease_date\tsha256")
        for row in rows:
            print(
                f"{row['genome_id']}\t{row['label']}\t{row['length_nt']}\t"
                f"{row['contig_count']}\t{row['release_date']}\t{row['sha256']}"
            )
    return 0


def _run_schema(args: argparse.Namespace) -> int:
    from .output import REFERENCE_COLUMNS, TRUTH_COLUMNS
    from .references import METADATA_RECOMMENDED_COLUMNS, METADATA_REQUIRED_COLUMNS
    from .schema_resources import load_schema

    header_contracts = {
        "metadata": (*METADATA_REQUIRED_COLUMNS, *METADATA_RECOMMENDED_COLUMNS),
        "truth": TRUTH_COLUMNS,
        "references": REFERENCE_COLUMNS,
    }
    columns = header_contracts.get(args.name)
    if columns is not None:
        print("\t".join(columns))
        return 0
    print(json.dumps(load_schema(args.name), ensure_ascii=False, indent=2, sort_keys=True))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Run CHIMERA and return a documented process status code."""

    parser = build_parser()
    args = parser.parse_args(argv)
    level = logging.ERROR if args.quiet else getattr(logging, args.log_level)
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s", stream=sys.stderr)
    try:
        if args.command in {"suite", "generate"}:
            return _run_generate(args)
        if args.command == "validate":
            return _run_validate(args)
        if args.command == "inspect":
            return _run_inspect(args)
        if args.command == "schema":
            return _run_schema(args)
    except IntegrityError as exc:
        LOGGER.error("Integrity validation failed: %s", exc)
        return 3
    except (ChimeraError, ValueError, OSError) as exc:
        LOGGER.error("%s", exc)
        return 2
    parser.error(f"unsupported command {args.command!r}")


def legacy_main() -> int:
    """Compatibility entry point retaining the legacy executable name."""

    print(
        "WARNING: 'metagenome-generator' is deprecated; use 'chimera' instead.",
        file=sys.stderr,
    )
    return main()
