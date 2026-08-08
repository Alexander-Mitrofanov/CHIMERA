"""One-command orchestration of reproducible CHIMERA benchmark bundles."""

from __future__ import annotations

import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import UTC, date, datetime
from enum import Enum
from pathlib import Path

from . import __version__
from .config import BenchmarkConfig, SimilarityBands, SplitKind
from .errors import ConfigurationError, IntegrityError
from .fasta import discover_fasta_files, write_fasta
from .fragments import derive_seed, generate_fragments, split_fragments_random
from .models import Contig, Fragment, Genome, Label
from .output import (
    REFERENCE_COLUMNS,
    SEQUENCE_COLUMNS,
    TRUTH_COLUMNS,
    fragment_statistics,
    reference_rows,
    sequence_rows,
    truth_rows,
    write_checksums,
    write_json,
    write_text,
    write_tsv,
)
from .provenance import software_content_sha256
from .references import ReferenceCatalog, load_reference_catalog
from .schema_resources import JSON_SCHEMA_NAMES, load_schema, schema_filename
from .similarity import format_similarity_value
from .splits import (
    SplitAssignment,
    SplitPartition,
    SplitPlan,
    build_split_plan,
)
from .validation import validate_bundle

BUNDLE_SCHEMA = "urn:chimera:benchmark-bundle:2"
SPLIT_SCHEMA = "urn:chimera:split-manifest:2"
ASSIGNMENT_COLUMNS = (
    "genome_id",
    "group_id",
    "label",
    "partition",
    "candidate_partition",
    "reason",
    "release_date",
    "taxon",
    "similarity_bin",
    "nearest_train_genome_id",
    "max_train_similarity",
    "similarity_coverage",
    "similarity_method",
    "strict_gate_train_genome_id",
    "strict_gate_similarity",
    "strict_gate_coverage",
    "strict_gate_method",
)
EXCLUSION_COLUMNS = (
    "genome_id",
    "label",
    "split",
    "reason",
    "duplicate_of",
    "source_sha256",
    "source_accession_version",
    "release_date",
    "nearest_train_genome_id",
    "max_train_similarity",
    "similarity_coverage",
    "similarity_method",
    "strict_gate_train_genome_id",
    "strict_gate_similarity",
    "strict_gate_coverage",
    "strict_gate_method",
)


@dataclass(frozen=True, slots=True)
class BenchmarkResult:
    """Small machine-facing result returned by library and CLI entry points."""

    output_dir: Path
    summary: Mapping[str, object]
    dry_run: bool = False

    def as_dict(self) -> dict[str, object]:
        return {
            "dry_run": self.dry_run,
            "output_dir": str(self.output_dir),
            "summary": dict(self.summary),
        }


@dataclass(frozen=True, slots=True)
class _PreparedBenchmark:
    catalog: ReferenceCatalog
    plans: Mapping[SplitKind, SplitPlan]


def _jsonable(value: object) -> object:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (date, datetime)):
        return value.isoformat()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, SimilarityBands):
        return asdict(value)
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in sorted(value.items())}
    if isinstance(value, (tuple, list, set, frozenset)):
        return [_jsonable(item) for item in value]
    raise TypeError(f"Cannot serialize provenance value of type {type(value).__name__}")


def _safe_target(target: Path) -> Path:
    resolved = target.expanduser().resolve()
    home = Path.home().resolve()
    working_directory = Path.cwd().resolve()
    protected = {
        home,
        working_directory,
        *home.parents,
        *working_directory.parents,
    }
    if resolved in protected:
        raise ConfigurationError(
            f"Refusing to use broad/protected directory as --outdir: {resolved}"
        )
    if not resolved.name:
        raise ConfigurationError("--outdir must name a dedicated benchmark directory")
    return resolved


_PRIVATE_DIRECTORY_MODE = 0o700
_PUBLISHED_DIRECTORY_MODE = 0o755
_PUBLISHED_FILE_MODE = 0o644


def _ensure_public_directory(path: Path) -> None:
    """Create only missing directory components and publish them as ``0755``."""

    missing: list[Path] = []
    cursor = path
    while not cursor.exists():
        missing.append(cursor)
        parent = cursor.parent
        if parent == cursor:
            break
        cursor = parent
    if cursor.exists() and not cursor.is_dir():
        raise ConfigurationError(f"Output parent is not a directory: {cursor}")
    for directory in reversed(missing):
        try:
            directory.mkdir(mode=_PRIVATE_DIRECTORY_MODE)
        except FileExistsError as error:
            if not directory.is_dir():
                raise ConfigurationError(
                    f"Output parent is not a directory: {directory}"
                ) from error
        else:
            directory.chmod(_PUBLISHED_DIRECTORY_MODE)


def _publish_bundle_permissions(staging: Path) -> None:
    """Set exact shared-HPC modes while keeping the staging root private until last."""

    for path in sorted(staging.rglob("*"), key=lambda item: item.as_posix()):
        if path.is_symlink():
            raise IntegrityError(f"Refusing to publish symbolic link in bundle: {path}")
        if path.is_dir():
            path.chmod(_PUBLISHED_DIRECTORY_MODE)
        elif path.is_file():
            path.chmod(_PUBLISHED_FILE_MODE)
    staging.chmod(_PUBLISHED_DIRECTORY_MODE)


def _require_valid_existing_bundle(target: Path) -> None:
    """Require an exact marker and a fully valid bundle before replacement."""

    marker = target / ".chimera-bundle"
    expected_marker = f"{BUNDLE_SCHEMA}\n"
    try:
        marker_content = marker.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise ConfigurationError(
            f"Refusing to replace unrecognized directory {target}; "
            "its .chimera-bundle marker is missing or unreadable"
        ) from exc
    if marker_content != expected_marker:
        raise ConfigurationError(
            f"Refusing to replace unrecognized directory {target}; "
            f".chimera-bundle must contain exactly {BUNDLE_SCHEMA!r}"
        )
    try:
        validate_bundle(target)
    except (IntegrityError, OSError) as exc:
        raise ConfigurationError(
            f"Refusing to replace invalid CHIMERA bundle {target}; "
            f"run 'chimera validate {target}' and repair or move it aside: {exc}"
        ) from exc


@contextmanager
def _atomic_bundle_directory(target: Path, *, overwrite: bool) -> Iterator[Path]:
    """Build beside the destination, then atomically swap a recognized bundle."""

    target = _safe_target(target)
    _ensure_public_directory(target.parent)
    if target.exists():
        if not overwrite:
            raise FileExistsError(f"Output directory already exists: {target}; use --force")
        _require_valid_existing_bundle(target)
    staging = Path(tempfile.mkdtemp(prefix=f".{target.name}.staging-", dir=target.parent))
    staging.chmod(_PRIVATE_DIRECTORY_MODE)
    backup: Path | None = None
    committed = False
    try:
        write_text(staging / ".chimera-bundle", f"{BUNDLE_SCHEMA}\n")
        yield staging
        if target.exists():
            _require_valid_existing_bundle(target)
            backup = target.with_name(f".{target.name}.backup-{os.getpid()}")
            if backup.exists():
                raise ConfigurationError(f"Atomic backup path already exists: {backup}")
        _publish_bundle_permissions(staging)
        if backup is not None:
            target.replace(backup)
        staging.replace(target)
        committed = True
        if backup is not None:
            shutil.rmtree(backup)
    except BaseException:
        if backup is not None and backup.exists() and not target.exists():
            backup.replace(target)
        raise
    finally:
        if not committed and staging.exists():
            shutil.rmtree(staging)


def _prepare(config: BenchmarkConfig) -> _PreparedBenchmark:
    catalog = load_reference_catalog(
        config.virus_paths,
        config.host_paths,
        metadata_path=config.metadata_path,
        duplicate_policy=config.duplicate_policy,
    )
    counts = Counter(genome.label for genome in catalog.genomes)
    for label in Label:
        if counts[label] < 2:
            raise ConfigurationError(
                f"At least two independent {label.value} genomes are required; "
                f"found {counts[label]} after integrity filtering"
            )
    longest_unsupported = [
        (genome.genome_id, max(contig.length for contig in genome.contigs), length)
        for genome in catalog.genomes
        for length in config.fragment_lengths
        if max(contig.length for contig in genome.contigs) < length
    ]
    if longest_unsupported:
        genome_id, longest, requested = longest_unsupported[0]
        raise ConfigurationError(
            f"Genome {genome_id!r} has longest contig {longest} nt and cannot emit "
            f"{requested}-nt fragments; choose a shorter --fragment-length"
        )
    plans: dict[SplitKind, SplitPlan] = {}
    for kind in config.splits:
        if kind is SplitKind.RANDOM:
            continue
        plans[kind] = build_split_plan(
            kind,
            catalog.genomes,
            test_fraction=config.test_fraction,
            seed=derive_seed(config.seed, "split-plan", kind.value),
            missing_metadata=config.missing_metadata,
            temporal_cutoff=config.temporal_cutoff,
            taxonomy_rank=config.taxonomy_rank,
            holdout_taxa=config.holdout_taxa,
            auto_holdout_count=config.auto_holdout_count,
            similarity_k=config.similarity_k,
            sketch_size=config.sketch_size,
            max_train_similarity=config.max_train_similarity,
            min_similarity_coverage=config.min_similarity_coverage,
            similarity_bands=config.similarity_bands,
            similarity_table=config.similarity_table,
        )
    return _PreparedBenchmark(catalog=catalog, plans=plans)


def _semantic_order(
    fragments: Iterable[Fragment], *, seed: int, namespace: str
) -> tuple[Fragment, ...]:
    return tuple(
        sorted(
            fragments,
            key=lambda fragment: (
                derive_seed(seed, "output-order", namespace, fragment.fragment_id),
                fragment.fragment_id,
            ),
        )
    )


def _assignment_row(assignment: SplitAssignment) -> dict[str, object]:
    return {
        "genome_id": assignment.genome_id,
        "group_id": assignment.group_id or assignment.genome_id,
        "label": assignment.label.value,
        "partition": assignment.partition.value,
        "candidate_partition": (
            assignment.candidate_partition.value if assignment.candidate_partition else ""
        ),
        "reason": assignment.reason,
        "release_date": assignment.release_date.isoformat() if assignment.release_date else "",
        "taxon": assignment.taxon or "",
        "similarity_bin": assignment.similarity_bin or "",
        "nearest_train_genome_id": assignment.nearest_train_genome_id or "",
        "max_train_similarity": (
            ""
            if assignment.max_train_similarity is None
            else format_similarity_value(assignment.max_train_similarity)
        ),
        "similarity_coverage": (
            ""
            if assignment.similarity_coverage is None
            else format_similarity_value(assignment.similarity_coverage)
        ),
        "similarity_method": assignment.similarity_method or "",
        "strict_gate_train_genome_id": assignment.strict_gate_train_genome_id or "",
        "strict_gate_similarity": (
            ""
            if assignment.strict_gate_similarity is None
            else format_similarity_value(assignment.strict_gate_similarity)
        ),
        "strict_gate_coverage": (
            ""
            if assignment.strict_gate_coverage is None
            else format_similarity_value(assignment.strict_gate_coverage)
        ),
        "strict_gate_method": assignment.strict_gate_method or "",
    }


def _assignment_details(plan: SplitPlan) -> dict[str, dict[str, object]]:
    return {
        assignment.genome_id: {
            "similarity_bin": assignment.similarity_bin or "",
            "max_train_similarity": assignment.max_train_similarity,
            "nearest_train_genome_id": assignment.nearest_train_genome_id or "",
        }
        for assignment in plan.assignments
    }


def _random_assignment_rows(genomes: Iterable[Genome]) -> list[dict[str, object]]:
    return [
        {
            "genome_id": genome.genome_id,
            "group_id": f"sha256:{genome.digest}",
            "label": genome.label.value,
            "partition": "both",
            "candidate_partition": "",
            "reason": "diagnostic_random_fragment_split",
            "release_date": (
                genome.metadata.release_date.isoformat() if genome.metadata.release_date else ""
            ),
            "taxon": "",
            "similarity_bin": "",
            "nearest_train_genome_id": "",
            "max_train_similarity": "",
            "similarity_coverage": "",
            "similarity_method": "",
            "strict_gate_train_genome_id": "",
            "strict_gate_similarity": "",
            "strict_gate_coverage": "",
            "strict_gate_method": "",
        }
        for genome in sorted(genomes, key=lambda item: item.genome_id)
    ]


def _split_validation(
    kind: SplitKind,
    train: tuple[Fragment, ...],
    test: tuple[Fragment, ...],
    genomes: Mapping[str, Genome],
    plan: SplitPlan | None,
) -> dict[str, object]:
    train_ids = {fragment.fragment_id for fragment in train}
    test_ids = {fragment.fragment_id for fragment in test}
    if train_ids & test_ids:
        raise IntegrityError(f"{kind.value}: fragment identifiers overlap train/test")
    train_groups = {fragment.genome_id for fragment in train}
    test_groups = {fragment.genome_id for fragment in test}
    shared_groups = train_groups & test_groups
    if kind is SplitKind.RANDOM:
        if not shared_groups:
            raise IntegrityError("Random diagnostic unexpectedly has no source-genome overlap")
    elif shared_groups:
        raise IntegrityError(
            f"{kind.value}: source genomes leak across train/test: {sorted(shared_groups)[:5]}"
        )
    for label in Label:
        if not any(fragment.label is label for fragment in train):
            raise IntegrityError(f"{kind.value}: training partition lacks {label.value}")
        if not any(fragment.label is label for fragment in test):
            raise IntegrityError(f"{kind.value}: test partition lacks {label.value}")
    shared_hashes = {genomes[group].digest for group in train_groups} & {
        genomes[group].digest for group in test_groups
    }
    if kind is not SplitKind.RANDOM and shared_hashes:
        raise IntegrityError(f"{kind.value}: canonical genome content leaks across partitions")
    if plan is not None:
        plan.validate()
    train_fragment_hashes = {fragment.digest for fragment in train}
    test_fragment_hashes = {fragment.digest for fragment in test}
    shared_fragment_hashes = train_fragment_hashes & test_fragment_hashes
    if kind is SplitKind.SIMILARITY and shared_fragment_hashes:
        raise IntegrityError(
            "similarity: exact fragment content remains in both strict train and test"
        )
    coordinate_overlap = _test_fragments_with_coordinate_overlap(train, test, genomes)
    return {
        "status": "pass",
        "diagnostic_only": kind is SplitKind.RANDOM,
        "fragment_id_overlap": len(train_ids & test_ids),
        "exact_fragment_content_overlap": len(shared_fragment_hashes),
        "test_fragments_with_coordinate_overlap": coordinate_overlap,
        "source_genome_overlap": len(shared_groups),
        "source_content_hash_overlap": len(shared_hashes),
        "train_source_genomes": len(train_groups),
        "test_source_genomes": len(test_groups),
    }


def _test_fragments_with_coordinate_overlap(
    train: Iterable[Fragment],
    test: Iterable[Fragment],
    genomes: Mapping[str, Genome],
) -> int:
    """Count test fragments whose source interval overlaps any train interval."""

    source_layout = {
        (genome.genome_id, contig.sequence_id): (contig.length, contig.topology)
        for genome in genomes.values()
        for contig in genome.contigs
    }

    def segments(fragment: Fragment) -> tuple[tuple[int, int], ...]:
        length, topology = source_layout[(fragment.genome_id, fragment.sequence_id)]
        if topology == "linear" or fragment.end <= length:
            return ((fragment.start, fragment.end),)
        return ((fragment.start, length), (0, fragment.end - length))

    intervals: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for fragment in train:
        intervals.setdefault((fragment.genome_id, fragment.sequence_id), []).extend(
            segments(fragment)
        )
    merged_by_source: dict[tuple[str, str], tuple[tuple[int, int], ...]] = {}
    for source, values in intervals.items():
        merged: list[list[int]] = []
        for start, end in sorted(values):
            if not merged or start >= merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        merged_by_source[source] = tuple((start, end) for start, end in merged)
    count = 0
    for fragment in test:
        source_intervals = merged_by_source.get((fragment.genome_id, fragment.sequence_id), ())
        if any(
            train_start < test_end and test_start < train_end
            for train_start, train_end in source_intervals
            for test_start, test_end in segments(fragment)
        ):
            count += 1
    return count


def _write_partition(
    split_dir: Path,
    name: str,
    fragments: tuple[Fragment, ...],
    *,
    genomes: Mapping[str, Genome],
    details: Mapping[str, Mapping[str, object]],
    semantic_partitions: Mapping[str, str] | None = None,
    view: str | None = None,
) -> list[dict[str, object]]:
    write_fasta(fragments, split_dir / f"{name}.fasta.gz")
    rows = truth_rows(
        fragments,
        partition=semantic_partitions or name,
        view=view,
        genomes=genomes,
        assignment_details=details,
    )
    write_tsv(split_dir / f"{name}.truth.tsv.gz", rows, TRUTH_COLUMNS)
    return rows


def _write_split(
    root: Path,
    config: BenchmarkConfig,
    kind: SplitKind,
    prepared: _PreparedBenchmark,
    input_digests: Mapping[Path, str],
) -> dict[str, object]:
    split_dir = root / kind.directory_name
    split_dir.mkdir(parents=True)
    genomes_by_id = prepared.catalog.by_id()
    plan = prepared.plans.get(kind)
    details: dict[str, dict[str, object]] = {}

    if kind is SplitKind.RANDOM:
        assignment_rows = _random_assignment_rows(prepared.catalog.genomes)
        write_tsv(split_dir / "assignments.tsv", assignment_rows, ASSIGNMENT_COLUMNS)
        all_fragments = generate_fragments(
            prepared.catalog.genomes,
            fragment_lengths=config.fragment_lengths,
            fragments_per_genome=config.fragments_per_genome,
            seed=derive_seed(config.seed, "fragments", "shared-library-v1"),
            strand_mode=config.strand_mode,
            max_ambiguous_fraction=config.max_ambiguous_fraction,
        )
        train, test = split_fragments_random(
            all_fragments,
            test_fraction=config.test_fraction,
            seed=derive_seed(config.seed, "fragment-membership", kind.value),
        )
        candidate_test: tuple[Fragment, ...] = ()
    else:
        assert plan is not None
        assignment_rows = [_assignment_row(assignment) for assignment in plan.assignments]
        # The immutable source plan is persisted before sequence generation.
        write_tsv(split_dir / "assignments.tsv", assignment_rows, ASSIGNMENT_COLUMNS)
        details = _assignment_details(plan)
        active_ids = plan.train_ids | plan.test_ids
        if kind is SplitKind.SIMILARITY:
            active_ids |= frozenset(
                assignment.genome_id
                for assignment in plan.excluded
                if assignment.candidate_partition is SplitPartition.TEST
            )
        active = tuple(genomes_by_id[genome_id] for genome_id in sorted(active_ids))
        all_fragments = generate_fragments(
            active,
            fragment_lengths=config.fragment_lengths,
            fragments_per_genome=config.fragments_per_genome,
            seed=derive_seed(config.seed, "fragments", "shared-library-v1"),
            strand_mode=config.strand_mode,
            max_ambiguous_fraction=config.max_ambiguous_fraction,
        )
        train = tuple(
            fragment for fragment in all_fragments if fragment.genome_id in plan.train_ids
        )
        test = tuple(fragment for fragment in all_fragments if fragment.genome_id in plan.test_ids)
        candidate_ids = {
            assignment.genome_id
            for assignment in plan.assignments
            if assignment.candidate_partition is SplitPartition.TEST
        }
        candidate_test = tuple(
            fragment for fragment in all_fragments if fragment.genome_id in candidate_ids
        )

    train = _semantic_order(train, seed=config.seed, namespace=f"{kind.value}:train")
    test = _semantic_order(test, seed=config.seed, namespace=f"{kind.value}:test")
    validation = _split_validation(kind, train, test, genomes_by_id, plan)
    train_truth = _write_partition(
        split_dir, "train", train, genomes=genomes_by_id, details=details
    )
    test_truth = _write_partition(split_dir, "test", test, genomes=genomes_by_id, details=details)

    exclusions: list[dict[str, object]] = []
    if plan is not None:
        exclusions.extend(
            {
                "genome_id": assignment.genome_id,
                "label": assignment.label.value,
                "split": kind.value,
                "reason": assignment.reason,
                "duplicate_of": "",
                "source_sha256": genomes_by_id[assignment.genome_id].digest,
                "source_accession_version": (
                    genomes_by_id[assignment.genome_id].metadata.accession_version or ""
                ),
                "release_date": (
                    assignment.release_date.isoformat()
                    if assignment.release_date is not None
                    else ""
                ),
                "nearest_train_genome_id": assignment.nearest_train_genome_id or "",
                "max_train_similarity": (
                    ""
                    if assignment.max_train_similarity is None
                    else format_similarity_value(assignment.max_train_similarity)
                ),
                "similarity_coverage": (
                    ""
                    if assignment.similarity_coverage is None
                    else format_similarity_value(assignment.similarity_coverage)
                ),
                "similarity_method": assignment.similarity_method or "",
                "strict_gate_train_genome_id": assignment.strict_gate_train_genome_id or "",
                "strict_gate_similarity": (
                    ""
                    if assignment.strict_gate_similarity is None
                    else format_similarity_value(assignment.strict_gate_similarity)
                ),
                "strict_gate_coverage": (
                    ""
                    if assignment.strict_gate_coverage is None
                    else format_similarity_value(assignment.strict_gate_coverage)
                ),
                "strict_gate_method": assignment.strict_gate_method or "",
            }
            for assignment in plan.excluded
        )
    write_tsv(split_dir / "excluded.tsv", exclusions, EXCLUSION_COLUMNS)

    if kind is SplitKind.SIMILARITY:
        assert plan is not None
        if config.similarity_table is not None:
            try:
                external_similarity_text = config.similarity_table.read_text(encoding="utf-8")
            except (OSError, UnicodeError) as error:
                raise ConfigurationError(
                    f"Cannot snapshot external similarity table {config.similarity_table}: {error}"
                ) from error
            write_text(split_dir / "external-similarity.tsv", external_similarity_text)
        candidate_partitions = {
            assignment.genome_id: assignment.partition.value
            for assignment in plan.assignments
            if assignment.candidate_partition is SplitPartition.TEST
        }
        candidate_test = _semantic_order(
            candidate_test, seed=config.seed, namespace=f"{kind.value}:candidate-test"
        )
        _write_partition(
            split_dir,
            "candidate_test",
            candidate_test,
            genomes=genomes_by_id,
            details=details,
            semantic_partitions=candidate_partitions,
            view="candidate_test",
        )
        strata_dir = split_dir / "test_strata"
        strata_dir.mkdir()
        for similarity_bin in (
            "high_similarity",
            "moderate_similarity",
            "low_similarity",
            "distant_detectable",
            "no_detectable_match",
        ):
            stratum = tuple(
                fragment
                for fragment in candidate_test
                if details[fragment.genome_id]["similarity_bin"] == similarity_bin
            )
            _write_partition(
                strata_dir,
                similarity_bin,
                stratum,
                genomes=genomes_by_id,
                details=details,
                semantic_partitions=candidate_partitions,
                view=f"test_strata/{similarity_bin}",
            )

    parameters: Mapping[str, object]
    if plan is None:
        parameters = {
            "diagnostic_only": True,
            "split_unit": "fragment",
            "test_fraction": config.test_fraction,
        }
    else:
        parameters = plan.parameters
    if kind is SplitKind.SIMILARITY and config.similarity_table is not None:
        similarity_source = _content_identifier(input_digests[config.similarity_table.resolve()])
        portable_parameters = dict(parameters)
        portable_parameters["similarity_source"] = similarity_source
        portable_parameters["similarity_table"] = similarity_source
        parameters = portable_parameters
    manifest = {
        "schema": SPLIT_SCHEMA,
        "protocol": kind.value,
        "protocol_id": kind.value,
        "parameters": _jsonable(parameters),
        "validation": validation,
        "train": fragment_statistics(train),
        "test": fragment_statistics(test),
        "truth_rows": {"train": len(train_truth), "test": len(test_truth)},
        "excluded_genomes": len(exclusions),
    }
    if kind is SplitKind.SIMILARITY:
        manifest["candidate_test"] = fragment_statistics(candidate_test)
    write_json(split_dir / "split.json", manifest)
    return manifest


def _root_exclusion_rows(catalog: ReferenceCatalog) -> list[dict[str, object]]:
    return [
        {
            "genome_id": exclusion.genome_id,
            "label": exclusion.label.value,
            "split": "reference_preflight",
            "reason": exclusion.reason,
            "duplicate_of": exclusion.duplicate_of or "",
            "source_sha256": exclusion.source_sha256,
            "source_accession_version": exclusion.accession_version or "",
            "release_date": exclusion.release_date.isoformat() if exclusion.release_date else "",
            "nearest_train_genome_id": "",
            "max_train_similarity": "",
            "similarity_coverage": "",
            "similarity_method": "",
            "strict_gate_train_genome_id": "",
            "strict_gate_similarity": "",
            "strict_gate_coverage": "",
            "strict_gate_method": "",
        }
        for exclusion in catalog.exclusions
    ]


def _report_markdown(
    config: BenchmarkConfig,
    prepared: _PreparedBenchmark,
    split_manifests: Mapping[str, Mapping[str, object]],
) -> str:
    lines = [
        "# CHIMERA benchmark report",
        "",
        f"Generated with CHIMERA {__version__} and master seed `{config.seed}`.",
        "All fragment identifiers are opaque; labels and coordinates live only in truth tables.",
        "",
        "## Reference inventory",
        "",
    ]
    counts = Counter(genome.label.value for genome in prepared.catalog.genomes)
    lines.append(
        f"Validated source genomes: **{len(prepared.catalog.genomes)}** "
        f"({counts['virus']} virus, {counts['host']} host)."
    )
    lines.extend(["", "## Evaluation protocols", ""])
    descriptions = {
        SplitKind.RANDOM.value: "Diagnostic only: fragments from every source occur on both sides.",
        SplitKind.GENOME.value: "Complete source genomes are held out before fragmentation.",
        SplitKind.SIMILARITY.value: "Genome-disjoint candidates are stratified by train similarity; `test` is strict.",
        SplitKind.TEMPORAL.value: "Training uses first public release dates on/before the inclusive cutoff.",
        SplitKind.TAXONOMY.value: "Selected viral taxa are absent from training; host negatives remain genome-disjoint.",
    }
    for name, manifest in split_manifests.items():
        stats = manifest.get("test")
        if not isinstance(stats, Mapping):
            raise IntegrityError(f"Split manifest {name!r} lacks test statistics")
        lines.append(
            f"- `{name}` — {descriptions[name]} Test records: {stats['records']}; "
            f"source genomes: {stats['source_genomes']}."
        )
    lines.extend(
        [
            "",
            "## Interpretation limits",
            "",
            "- Random-fragment performance does not measure unseen-genome generalization.",
            "- Built-in similarity values are Mash-style MinHash estimates, not universal viral taxonomy boundaries.",
            "- A temporal run without an archived historical snapshot is retrospective release-date filtering.",
            "- Taxonomic novelty and sequence novelty are distinct; report both where possible.",
            "- These are exact synthetic DNA fragments, not empirically calibrated sequencing reads.",
        ]
    )
    return "\n".join(lines)


def generate_benchmark(config: BenchmarkConfig, *, dry_run: bool = False) -> BenchmarkResult:
    """Preflight, plan, generate, validate, and atomically publish a benchmark suite."""

    if not isinstance(config, BenchmarkConfig):
        raise TypeError("config must be a BenchmarkConfig")
    output_dir = _safe_target(config.output_dir)
    if output_dir.exists() and not config.overwrite and not dry_run:
        raise FileExistsError(f"Output directory already exists: {output_dir}; use --force")
    prepared = _prepare(config)
    preflight_summary: dict[str, object] = {
        "genomes": len(prepared.catalog.genomes),
        "host_genomes": sum(genome.label is Label.HOST for genome in prepared.catalog.genomes),
        "protocols": [kind.value for kind in config.splits],
        "virus_genomes": sum(genome.label is Label.VIRUS for genome in prepared.catalog.genomes),
    }
    if dry_run:
        preflight_summary["resolved_plans"] = {
            kind.value: {
                "train": len(plan.train),
                "test": len(plan.test),
                "excluded": len(plan.excluded),
                "parameters": _jsonable(plan.parameters),
            }
            for kind, plan in prepared.plans.items()
        }
        return BenchmarkResult(output_dir=output_dir, summary=preflight_summary, dry_run=True)

    started = datetime.now(UTC)
    input_digests = _input_digest_map(config, prepared.catalog)
    with _atomic_bundle_directory(output_dir, overwrite=config.overwrite) as staging:
        resolved_config = _semantic_config(config, input_digests)
        write_json(staging / "resolved-config.json", resolved_config)
        schema_dir = staging / "schemas"
        schema_dir.mkdir()
        for schema_name in JSON_SCHEMA_NAMES:
            write_json(
                schema_dir / schema_filename(schema_name),
                load_schema(schema_name),
            )
        references = reference_rows(
            prepared.catalog.genomes,
            source_digests=input_digests,
        )
        write_tsv(staging / "references.tsv", references, REFERENCE_COLUMNS)
        sequences = sequence_rows(
            prepared.catalog.genomes,
            source_digests=input_digests,
        )
        write_tsv(staging / "sequences.tsv", sequences, SEQUENCE_COLUMNS)
        source_contigs = tuple(
            Contig(
                sequence_id=contig.sequence_id,
                sequence=contig.sequence,
                accession_version=contig.accession_version,
                release_date=contig.release_date,
                topology=contig.topology,
                taxonomy=contig.taxonomy,
                metadata_extra=contig.metadata_extra,
            )
            for genome in sorted(prepared.catalog.genomes, key=lambda item: item.genome_id)
            for contig in sorted(genome.contigs, key=lambda item: item.sequence_id)
        )
        write_fasta(source_contigs, staging / "source-sequences.fasta.gz")
        write_tsv(
            staging / "excluded.tsv",
            _root_exclusion_rows(prepared.catalog),
            EXCLUSION_COLUMNS,
        )
        split_manifests: dict[str, Mapping[str, object]] = {}
        for kind in config.splits:
            split_manifests[kind.value] = _write_split(
                staging,
                config,
                kind,
                prepared,
                input_digests,
            )
        manifest = {
            "schema": BUNDLE_SCHEMA,
            "tool": {
                "name": "CHIMERA",
                "version": __version__,
                "software_content_sha256": software_content_sha256(),
                **_git_provenance(),
            },
            "data_model": {
                "alphabet": "IUPAC DNA",
                "coordinate_system": "0-based-half-open",
                "coordinate_systems": {
                    "linear": "0-based-half-open",
                    "circular": "0-based-half-open-circular",
                },
                "coordinate_semantics": {
                    "linear": (
                        "source_start is inclusive and source_end is exclusive; "
                        "0 <= source_start < source_end <= source_length"
                    ),
                    "circular": (
                        "coordinates are an unwrapped forward-source interval; source_start "
                        "is in [0, source_length), source_end equals source_start plus "
                        "fragment_length, and source_end may exceed source_length when the "
                        "interval wraps the origin"
                    ),
                },
                "fragment_headers": "opaque label-free identifiers",
                "grouping": "canonical_topology_aware_genome_sha256_v2",
                "synthetic": True,
            },
            "randomness": {
                "master_seed": config.seed,
                "algorithm": "Python random.Random with semantic BLAKE2b-derived sub-seeds",
                "seed_derivation": "chimera.seed.v1",
                "python_implementation": platform.python_implementation(),
                "python_version": platform.python_version(),
            },
            "references": {
                "count": len(prepared.catalog.genomes),
                "inputs": _input_inventory(config, prepared.catalog, input_digests),
                "preflight_exclusions": len(prepared.catalog.exclusions),
            },
            "splits": split_manifests,
        }
        write_json(staging / "manifest.json", manifest)
        write_text(staging / "REPORT.md", _report_markdown(config, prepared, split_manifests))
        write_json(
            staging / "execution.json",
            {
                "started_at_utc": started.isoformat(),
                "finished_at_utc": datetime.now(UTC).isoformat(),
                "python": sys.version.split()[0],
                "platform": platform.platform(),
                "command": [Path(sys.argv[0]).name, "<arguments recorded in resolved-config.json>"],
                "status": "complete",
            },
        )
        write_checksums(staging)
        validate_bundle(staging)
    return BenchmarkResult(output_dir=output_dir, summary=preflight_summary)


def _file_sha256(path: Path) -> str:
    from .output import sha256_file

    return sha256_file(path)


def _content_identifier(digest: str) -> str:
    return f"sha256:{digest}"


def _input_digest_map(
    config: BenchmarkConfig,
    catalog: ReferenceCatalog,
) -> dict[Path, str]:
    paths = set(catalog.source_files)
    if config.metadata_path is not None:
        paths.add(config.metadata_path.resolve())
    if config.similarity_table is not None:
        paths.add(config.similarity_table.resolve())
    return {
        path.resolve(): _file_sha256(path.resolve())
        for path in sorted(paths, key=lambda item: item.as_posix())
    }


def _semantic_config(
    config: BenchmarkConfig,
    input_digests: Mapping[Path, str],
) -> dict[str, object]:
    values = config.as_manifest_dict()
    values.pop("overwrite", None)
    values["output_dir"] = "bundle"

    def identifiers(paths: tuple[Path, ...]) -> list[str]:
        files = discover_fasta_files(paths)
        return sorted(_content_identifier(input_digests[path.resolve()]) for path in files)

    values["virus_paths"] = identifiers(config.virus_paths)
    values["host_paths"] = identifiers(config.host_paths)
    values["metadata_path"] = (
        _content_identifier(input_digests[config.metadata_path.resolve()])
        if config.metadata_path is not None
        else None
    )
    values["similarity_table"] = (
        _content_identifier(input_digests[config.similarity_table.resolve()])
        if config.similarity_table is not None
        else None
    )
    return values


def _input_inventory(
    config: BenchmarkConfig,
    catalog: ReferenceCatalog,
    input_digests: Mapping[Path, str],
) -> list[dict[str, str]]:
    roles: dict[Path, str] = dict.fromkeys(catalog.source_files, "reference_fasta")
    if config.metadata_path is not None:
        roles[config.metadata_path.resolve()] = "reference_metadata"
    if config.similarity_table is not None:
        roles[config.similarity_table.resolve()] = "external_similarity_table"
    return [
        {
            "content_id": _content_identifier(input_digests[path.resolve()]),
            "role": roles[path],
            "sha256": input_digests[path.resolve()],
        }
        for path in sorted(
            roles,
            key=lambda item: (roles[item], input_digests[item.resolve()]),
        )
    ]


def _git_provenance() -> dict[str, object]:
    source_root = Path(__file__).resolve().parents[2]
    if not (source_root / ".git").exists():
        return {"git_revision": "unknown", "git_dirty": None}
    try:
        revision_result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=source_root,
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
        status_result = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=normal"],
            cwd=source_root,
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (FileNotFoundError, subprocess.SubprocessError):
        return {"git_revision": "unknown", "git_dirty": None}
    revision = revision_result.stdout.strip()
    if re.fullmatch(r"[0-9a-f]{40}", revision) is None:
        return {"git_revision": "unknown", "git_dirty": None}
    return {"git_revision": revision, "git_dirty": bool(status_result.stdout.strip())}


__all__ = ["BenchmarkResult", "generate_benchmark"]
