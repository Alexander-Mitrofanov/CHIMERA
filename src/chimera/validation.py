"""Independent validation of published CHIMERA benchmark bundles.

This module deliberately re-computes integrity facts from FASTA, truth, and
assignment files.  Recorded ``validation.status`` values are provenance only
and are never trusted as evidence that a bundle is sound.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter
from collections.abc import Mapping
from contextlib import suppress
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Any, NoReturn, TextIO, cast

from . import __version__
from .config import BenchmarkConfig, SimilarityBands, SplitKind
from .errors import ChimeraError, IntegrityError
from .fasta import FastaFormatError, read_fasta
from .fragments import derive_seed, generate_fragments, split_fragments_random
from .models import (
    Contig,
    Fragment,
    Genome,
    GenomeMetadata,
    Label,
    Topology,
    canonical_sequence_hash,
    deterministic_genome_hash,
    reverse_complement,
)
from .output import REFERENCE_COLUMNS, SEQUENCE_COLUMNS, TRUTH_COLUMNS, verify_checksums
from .provenance import software_content_sha256
from .schema_resources import (
    JSON_SCHEMA_NAMES,
    SchemaName,
    load_schema,
    schema_filename,
    validate_instance,
)
from .similarity import format_similarity_value
from .splits import SplitAssignment, SplitPartition, SplitPlan, build_split_plan

_BUNDLE_SCHEMA = "urn:chimera:benchmark-bundle:2"
_SPLIT_SCHEMA = "urn:chimera:split-manifest:2"
_ROOT_FILES = (
    ".chimera-bundle",
    "REPORT.md",
    "checksums.sha256",
    "excluded.tsv",
    "execution.json",
    "manifest.json",
    "references.tsv",
    "resolved-config.json",
    "sequences.tsv",
    "source-sequences.fasta.gz",
)
_ASSIGNMENT_COLUMNS = (
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
_EXCLUSION_COLUMNS = (
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
_SIMILARITY_BINS = (
    "high_similarity",
    "moderate_similarity",
    "low_similarity",
    "distant_detectable",
    "no_detectable_match",
)
_PRIMARY_SPLIT_FILES = frozenset(
    {
        "assignments.tsv",
        "excluded.tsv",
        "split.json",
        "train.fasta.gz",
        "train.truth.tsv.gz",
        "test.fasta.gz",
        "test.truth.tsv.gz",
    }
)


@dataclass(frozen=True, slots=True)
class ValidationReport:
    """Structured result returned after an entire bundle passes validation."""

    root: Path
    checksums_verified: int
    split_counts: tuple[tuple[str, int, int], ...]
    primary_fasta_records_verified: int
    primary_truth_rows_verified: int
    auxiliary_fasta_records_verified: int
    auxiliary_truth_rows_verified: int
    assignment_rows_verified: int
    checks: tuple[str, ...]

    @property
    def status(self) -> str:
        """Machine-facing overall status."""

        return "pass"

    @property
    def fasta_records_verified(self) -> int:
        """Total records checked across primary partitions and auxiliary views."""

        return self.primary_fasta_records_verified + self.auxiliary_fasta_records_verified

    @property
    def truth_rows_verified(self) -> int:
        """Total truth rows checked across primary partitions and auxiliary views."""

        return self.primary_truth_rows_verified + self.auxiliary_truth_rows_verified

    def as_dict(self) -> dict[str, object]:
        """Return a stable JSON-compatible representation."""

        return {
            "status": self.status,
            "root": str(self.root),
            "checksums_verified": self.checksums_verified,
            "fasta_records_verified": self.fasta_records_verified,
            "truth_rows_verified": self.truth_rows_verified,
            "primary_fasta_records_verified": self.primary_fasta_records_verified,
            "primary_truth_rows_verified": self.primary_truth_rows_verified,
            "auxiliary_fasta_records_verified": self.auxiliary_fasta_records_verified,
            "auxiliary_truth_rows_verified": self.auxiliary_truth_rows_verified,
            "assignment_rows_verified": self.assignment_rows_verified,
            "splits": {
                kind: {"train_records": train, "test_records": test}
                for kind, train, test in self.split_counts
            },
            "checks": list(self.checks),
        }

    def summary(self) -> str:
        """Return a concise human-readable validation summary."""

        return (
            f"Validated CHIMERA bundle {self.root}: {len(self.split_counts)} split(s), "
            f"{self.fasta_records_verified} FASTA record(s) "
            f"({self.primary_fasta_records_verified} primary), "
            f"{self.truth_rows_verified} truth row(s) "
            f"({self.primary_truth_rows_verified} primary), and "
            f"{self.checksums_verified} checksum(s)."
        )

    def __str__(self) -> str:
        return self.summary()


@dataclass(frozen=True, slots=True)
class _Reference:
    genome_id: str
    label: str
    accession_version: str
    release_date: str
    digest: str
    sequence_ids: frozenset[str]
    source_input_ids: frozenset[str]
    taxonomy: Mapping[str, str]
    contig_count: int
    length_nt: int
    metadata_extra: Mapping[str, str]


@dataclass(frozen=True, slots=True)
class _SourceSequence:
    sequence_id: str
    genome_id: str
    label: str
    accession_version: str
    release_date: str
    topology: str
    sequence: str
    sha256: str
    canonical_sha256: str
    source_input_id: str
    taxonomy: Mapping[str, str]
    metadata_extra: Mapping[str, str]


@dataclass(frozen=True, slots=True)
class _Partition:
    name: str
    rows: tuple[dict[str, str], ...]
    sequences: Mapping[str, str]
    ordered_ids: tuple[str, ...]

    @property
    def ids(self) -> frozenset[str]:
        return frozenset(self.sequences)

    @property
    def source_genomes(self) -> frozenset[str]:
        return frozenset(row["source_genome_id"] for row in self.rows)


def _fail(message: str) -> NoReturn:
    raise IntegrityError(message)


def _require_files(root: Path, relative_paths: tuple[str, ...]) -> None:
    missing = [relative for relative in relative_paths if not (root / relative).is_file()]
    if missing:
        _fail("Bundle is missing required file(s): " + ", ".join(missing))


def _read_json(path: Path) -> dict[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                _fail(f"{path}: duplicate JSON key {key!r}")
            value[key] = item
        return value

    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=reject_duplicates)
    except IntegrityError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        _fail(f"Cannot read JSON object {path}: {error}")
    if not isinstance(value, dict):
        _fail(f"{path}: top-level JSON value must be an object")
    return cast(dict[str, Any], value)


def _open_text(path: Path) -> TextIO:
    try:
        if path.name.lower().endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8", newline="")
        return path.open("rt", encoding="utf-8", newline="")
    except OSError as error:
        _fail(f"Cannot open table {path}: {error}")


def _read_tsv(path: Path, required: tuple[str, ...]) -> tuple[dict[str, str], ...]:
    try:
        with _open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames or []
            if len(header) != len(set(header)):
                _fail(f"{path}: TSV header contains duplicate columns")
            missing = sorted(set(required) - set(header))
            if missing:
                _fail(f"{path}: missing required column(s): {', '.join(missing)}")
            rows: list[dict[str, str]] = []
            for line_number, row in enumerate(reader, start=2):
                if None in row:
                    _fail(f"{path}:{line_number}: row has more fields than the header")
                rows.append({key: (value or "").strip() for key, value in row.items()})
            return tuple(rows)
    except IntegrityError:
        raise
    except (OSError, UnicodeError, csv.Error) as error:
        _fail(f"Cannot read TSV {path}: {error}")


def _validate_logical_row(
    row: Mapping[str, str],
    schema: SchemaName,
    *,
    integer_fields: tuple[str, ...] = (),
) -> None:
    """Validate a TSV row after decoding fields typed as integers by its schema."""

    logical: dict[str, object] = dict(row)
    for field in integer_fields:
        with suppress(ValueError):
            logical[field] = int(row[field])
    validate_instance(logical, schema)


def _parse_int(value: str, *, path: Path, line_number: int, field: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise IntegrityError(
            f"{path}:{line_number}: {field} must be an integer, got {value!r}"
        ) from error
    if str(parsed) != value:
        _fail(f"{path}:{line_number}: {field} must use canonical integer syntax")
    return parsed


def _parse_fraction(
    value: str,
    *,
    path: Path,
    line_number: int,
    field: str,
    allow_empty: bool,
) -> float | None:
    if not value and allow_empty:
        return None
    try:
        parsed = float(value)
    except ValueError as error:
        raise IntegrityError(f"{path}:{line_number}: {field} must be a number on [0, 1]") from error
    if not math.isfinite(parsed) or not 0.0 <= parsed <= 1.0:
        _fail(f"{path}:{line_number}: {field} must be finite and on [0, 1]")
    return parsed


def _parse_date(value: str, *, path: Path, line_number: int, field: str) -> date:
    try:
        return date.fromisoformat(value)
    except ValueError as error:
        raise IntegrityError(
            f"{path}:{line_number}: {field} must be ISO YYYY-MM-DD, got {value!r}"
        ) from error


def _parse_json_string_map(
    value: str, *, path: Path, line_number: int, field: str
) -> dict[str, str]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, item in pairs:
            if key in result:
                _fail(f"{path}:{line_number}: duplicate {field} key {key!r}")
            result[key] = item
        return result

    try:
        parsed = json.loads(value, object_pairs_hook=reject_duplicates)
    except json.JSONDecodeError as error:
        raise IntegrityError(f"{path}:{line_number}: {field} must be a JSON object") from error
    if not isinstance(parsed, dict) or any(
        not isinstance(key, str) or not isinstance(item, str) for key, item in parsed.items()
    ):
        _fail(f"{path}:{line_number}: {field} must map strings to strings")
    return cast(dict[str, str], parsed)


def _parse_json_string_list(
    value: str, *, path: Path, line_number: int, field: str
) -> tuple[str, ...]:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as error:
        raise IntegrityError(f"{path}:{line_number}: {field} must be a JSON array") from error
    if (
        not isinstance(parsed, list)
        or not parsed
        or any(not isinstance(item, str) or not item for item in parsed)
        or len(parsed) != len(set(parsed))
    ):
        _fail(f"{path}:{line_number}: {field} must be a nonempty array of unique strings")
    return tuple(cast(list[str], parsed))


def _read_references(root: Path, manifest: Mapping[str, Any]) -> dict[str, _Reference]:
    path = root / "references.tsv"
    rows = _read_tsv(path, REFERENCE_COLUMNS)
    for row in rows:
        _validate_logical_row(
            row,
            SchemaName.REFERENCE_ROW,
            integer_fields=("contig_count", "length_nt"),
        )
    references: dict[str, _Reference] = {}
    seen_hashes: dict[str, str] = {}
    for line_number, row in enumerate(rows, start=2):
        genome_id = row["genome_id"]
        if not genome_id:
            _fail(f"{path}:{line_number}: genome_id must not be empty")
        if genome_id in references:
            _fail(f"{path}:{line_number}: duplicate genome_id {genome_id!r}")
        try:
            Label(row["label"])
        except ValueError as error:
            raise IntegrityError(
                f"{path}:{line_number}: label must be 'virus' or 'host'"
            ) from error
        digest = row["sha256"]
        if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
            _fail(f"{path}:{line_number}: invalid canonical genome SHA-256")
        if digest in seen_hashes:
            _fail(
                f"{path}:{line_number}: duplicate canonical genome content for "
                f"{genome_id!r} and {seen_hashes[digest]!r}"
            )
        seen_hashes[digest] = genome_id
        taxonomy = _parse_json_string_map(
            row["taxonomy"], path=path, line_number=line_number, field="taxonomy"
        )
        metadata_extra = _parse_json_string_map(
            row["metadata_extra"],
            path=path,
            line_number=line_number,
            field="metadata_extra",
        )
        sequence_ids = _parse_json_string_list(
            row["sequence_ids"], path=path, line_number=line_number, field="sequence_ids"
        )
        source_input_ids = _parse_json_string_list(
            row["source_input_ids"],
            path=path,
            line_number=line_number,
            field="source_input_ids",
        )
        if any(re.fullmatch(r"sha256:[0-9a-f]{64}", item) is None for item in source_input_ids):
            _fail(f"{path}:{line_number}: source_input_ids must be content-addressed")
        contig_count = _parse_int(
            row["contig_count"], path=path, line_number=line_number, field="contig_count"
        )
        length_nt = _parse_int(
            row["length_nt"], path=path, line_number=line_number, field="length_nt"
        )
        if contig_count < 1 or length_nt < 1:
            _fail(f"{path}:{line_number}: contig_count and length_nt must be positive")
        references[genome_id] = _Reference(
            genome_id=genome_id,
            label=row["label"],
            accession_version=row["accession_version"],
            release_date=row["release_date"],
            digest=digest,
            sequence_ids=frozenset(sequence_ids),
            source_input_ids=frozenset(source_input_ids),
            taxonomy=taxonomy,
            contig_count=contig_count,
            length_nt=length_nt,
            metadata_extra=metadata_extra,
        )
        if not references[genome_id].sequence_ids:
            _fail(f"{path}:{line_number}: sequence_ids must not be empty")
        if len(references[genome_id].sequence_ids) != contig_count:
            _fail(f"{path}:{line_number}: contig_count disagrees with sequence_ids")
    reference_manifest = manifest.get("references")
    if not isinstance(reference_manifest, Mapping):
        _fail(f"{root / 'manifest.json'}: references must be an object")
    if reference_manifest.get("count") != len(references):
        _fail("manifest.json reference count does not match references.tsv")
    return references


def _read_source_sequences(
    root: Path, references: Mapping[str, _Reference]
) -> dict[str, _SourceSequence]:
    table_path = root / "sequences.tsv"
    fasta_path = root / "source-sequences.fasta.gz"
    rows = _read_tsv(table_path, SEQUENCE_COLUMNS)
    for row in rows:
        _validate_logical_row(row, SchemaName.SEQUENCE_ROW, integer_fields=("length_nt",))
    try:
        contigs = read_fasta(fasta_path)
    except (FastaFormatError, OSError, ValueError) as error:
        raise IntegrityError(f"Cannot validate source FASTA {fasta_path}: {error}") from error
    fasta_ids = tuple(contig.sequence_id for contig in contigs)
    row_ids = tuple(row["sequence_id"] for row in rows)
    if fasta_ids != row_ids:
        _fail("source-sequences.fasta.gz and sequences.tsv must have identical ordered IDs")
    fasta_by_id = {contig.sequence_id: contig.sequence for contig in contigs}
    sources: dict[str, _SourceSequence] = {}
    by_genome: dict[str, list[Contig]] = {}
    for line_number, row in enumerate(rows, start=2):
        sequence_id = row["sequence_id"]
        genome_id = row["genome_id"]
        if not sequence_id or sequence_id in sources:
            _fail(f"{table_path}:{line_number}: empty/duplicate sequence_id {sequence_id!r}")
        reference = references.get(genome_id)
        if reference is None:
            _fail(f"{table_path}:{line_number}: unknown genome_id {genome_id!r}")
        if row["label"] != reference.label:
            _fail(f"{table_path}:{line_number}: label disagrees with references.tsv")
        if sequence_id not in reference.sequence_ids:
            _fail(f"{table_path}:{line_number}: sequence is absent from its reference group")
        topology = row["topology"]
        if topology not in {"linear", "circular"}:
            _fail(f"{table_path}:{line_number}: topology must be linear or circular")
        length = _parse_int(
            row["length_nt"], path=table_path, line_number=line_number, field="length_nt"
        )
        sequence = fasta_by_id[sequence_id]
        if length != len(sequence):
            _fail(f"{table_path}:{line_number}: length_nt disagrees with source FASTA")
        exact_digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        canonical_digest = canonical_sequence_hash(sequence, circular=topology == "circular")
        if row["sha256"] != exact_digest or row["canonical_sha256"] != canonical_digest:
            _fail(f"{table_path}:{line_number}: source sequence digest is incorrect")
        if row["release_date"]:
            _parse_date(
                row["release_date"],
                path=table_path,
                line_number=line_number,
                field="release_date",
            )
        metadata_extra = _parse_json_string_map(
            row["metadata_extra"],
            path=table_path,
            line_number=line_number,
            field="metadata_extra",
        )
        taxonomy = _parse_json_string_map(
            row["taxonomy"],
            path=table_path,
            line_number=line_number,
            field="taxonomy",
        )
        source = _SourceSequence(
            sequence_id=sequence_id,
            genome_id=genome_id,
            label=row["label"],
            accession_version=row["accession_version"],
            release_date=row["release_date"],
            topology=topology,
            sequence=sequence,
            sha256=exact_digest,
            canonical_sha256=canonical_digest,
            source_input_id=row["source_input_id"],
            taxonomy=taxonomy,
            metadata_extra=metadata_extra,
        )
        if re.fullmatch(r"sha256:[0-9a-f]{64}", source.source_input_id) is None:
            _fail(f"{table_path}:{line_number}: source_input_id must be content-addressed")
        sources[sequence_id] = source
        by_genome.setdefault(genome_id, []).append(
            Contig(
                sequence_id=sequence_id,
                sequence=sequence,
                accession_version=row["accession_version"] or None,
                release_date=(
                    date.fromisoformat(row["release_date"]) if row["release_date"] else None
                ),
                topology=cast(Topology, topology),
                taxonomy=tuple(sorted(taxonomy.items())),
                metadata_extra=tuple(sorted(metadata_extra.items())),
            )
        )
    expected_ids = frozenset(
        sequence_id for reference in references.values() for sequence_id in reference.sequence_ids
    )
    if frozenset(sources) != expected_ids:
        _fail("sequences.tsv must inventory every retained source sequence exactly once")
    for genome_id, reference in references.items():
        genome_contigs = tuple(by_genome.get(genome_id, ()))
        if len(genome_contigs) != reference.contig_count:
            _fail(f"Source contig count disagrees for genome {genome_id!r}")
        if sum(contig.length for contig in genome_contigs) != reference.length_nt:
            _fail(f"Source length disagrees for genome {genome_id!r}")
        if deterministic_genome_hash(genome_contigs) != reference.digest:
            _fail(f"Source content digest disagrees for genome {genome_id!r}")
        expected_inputs = {
            sources[sequence_id].source_input_id for sequence_id in reference.sequence_ids
        }
        if reference.source_input_ids != expected_inputs:
            _fail(f"Reference source input IDs disagree for genome {genome_id!r}")
        known_dates = tuple(contig.release_date for contig in genome_contigs)
        expected_release = (
            max(cast(tuple[date, ...], known_dates)).isoformat()
            if all(value is not None for value in known_dates)
            else ""
        )
        if reference.release_date != expected_release:
            _fail(
                f"Genome-level release_date for {genome_id!r} is not the latest "
                "fully known segment release date"
            )
        accessions = {
            contig.accession_version
            for contig in genome_contigs
            if contig.accession_version is not None
        }
        expected_accession = next(iter(accessions)) if len(accessions) == 1 else ""
        if reference.accession_version != expected_accession:
            _fail(
                f"Genome-level accession_version for {genome_id!r} disagrees with "
                "its per-sequence inventory"
            )
        taxonomy_by_rank: dict[str, set[str]] = {}
        for contig in genome_contigs:
            for rank, value in contig.taxonomy:
                taxonomy_by_rank.setdefault(rank, set()).add(value)
        if any(len(values) != 1 for values in taxonomy_by_rank.values()):
            _fail(f"Per-sequence taxonomy conflicts within genome {genome_id!r}")
        expected_taxonomy = {rank: next(iter(values)) for rank, values in taxonomy_by_rank.items()}
        if dict(reference.taxonomy) != expected_taxonomy:
            _fail(f"Genome taxonomy for {genome_id!r} disagrees with per-sequence inventory")
    return sources


def _reconstruct_genomes(
    references: Mapping[str, _Reference],
    sources: Mapping[str, _SourceSequence],
) -> tuple[Genome, ...]:
    """Reconstruct the exact retained genome catalog from checked bundle data."""

    genomes: list[Genome] = []
    for genome_id, reference in sorted(references.items()):
        contigs = tuple(
            Contig(
                sequence_id=source.sequence_id,
                sequence=source.sequence,
                accession_version=source.accession_version or None,
                release_date=(
                    date.fromisoformat(source.release_date) if source.release_date else None
                ),
                topology=cast(Topology, source.topology),
                taxonomy=tuple(sorted(source.taxonomy.items())),
                metadata_extra=tuple(sorted(source.metadata_extra.items())),
            )
            for source in sorted(sources.values(), key=lambda item: item.sequence_id)
            if source.genome_id == genome_id
        )
        genomes.append(
            Genome(
                genome_id=genome_id,
                label=Label(reference.label),
                contigs=contigs,
                metadata=GenomeMetadata(
                    release_date=(
                        date.fromisoformat(reference.release_date)
                        if reference.release_date
                        else None
                    ),
                    taxonomy=tuple(sorted(reference.taxonomy.items())),
                    accession_version=reference.accession_version or None,
                    extra=tuple(sorted(reference.metadata_extra.items())),
                ),
            )
        )
    return tuple(genomes)


def _read_assignments(
    path: Path,
    kind: SplitKind,
    references: Mapping[str, _Reference],
) -> dict[str, dict[str, str]]:
    rows = _read_tsv(path, _ASSIGNMENT_COLUMNS)
    for row in rows:
        _validate_logical_row(row, SchemaName.ASSIGNMENT_ROW)
    assignments: dict[str, dict[str, str]] = {}
    allowed_partitions = {"both"} if kind is SplitKind.RANDOM else {"train", "test", "excluded"}
    for line_number, row in enumerate(rows, start=2):
        genome_id = row["genome_id"]
        if genome_id in assignments:
            _fail(f"{path}:{line_number}: duplicate genome assignment {genome_id!r}")
        reference = references.get(genome_id)
        if reference is None:
            _fail(f"{path}:{line_number}: assignment names unknown genome {genome_id!r}")
        if row["partition"] not in allowed_partitions:
            _fail(
                f"{path}:{line_number}: partition {row['partition']!r} is invalid for "
                f"the {kind.value} protocol"
            )
        if row["label"] != reference.label:
            _fail(f"{path}:{line_number}: assignment label disagrees with references.tsv")
        expected_group = f"sha256:{reference.digest}"
        if row["group_id"] != expected_group:
            _fail(f"{path}:{line_number}: canonical content group disagrees with references.tsv")
        if row["release_date"] != reference.release_date:
            _fail(f"{path}:{line_number}: release_date disagrees with references.tsv")
        if row["similarity_bin"] and row["similarity_bin"] not in _SIMILARITY_BINS:
            _fail(f"{path}:{line_number}: unknown similarity_bin {row['similarity_bin']!r}")
        for field in (
            "max_train_similarity",
            "similarity_coverage",
            "strict_gate_similarity",
            "strict_gate_coverage",
        ):
            parsed = _parse_fraction(
                row[field],
                path=path,
                line_number=line_number,
                field=field,
                allow_empty=True,
            )
            if parsed is not None and row[field] != format_similarity_value(parsed):
                _fail(f"{path}:{line_number}: {field} must use canonical round-trip syntax")
        assignments[genome_id] = row
    if set(assignments) != set(references):
        missing = sorted(set(references) - set(assignments))
        extra = sorted(set(assignments) - set(references))
        _fail(
            f"{path}: assignments must cover every reference exactly once; "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )
    return assignments


def _read_partition(
    split_dir: Path,
    name: str,
    *,
    sources: Mapping[str, _SourceSequence],
    references: Mapping[str, _Reference],
    expected_view: str | None = None,
    allowed_partitions: frozenset[str] | None = None,
    allow_empty: bool = False,
) -> _Partition:
    """Read one physical FASTA/truth view and authenticate every source slice."""

    fasta_path = split_dir / f"{name}.fasta.gz"
    truth_path = split_dir / f"{name}.truth.tsv.gz"
    _require_files(split_dir, (f"{name}.fasta.gz", f"{name}.truth.tsv.gz"))
    try:
        contigs = read_fasta(fasta_path)
    except FastaFormatError as error:
        genuinely_empty = False
        if allow_empty:
            try:
                with _open_text(fasta_path) as handle:
                    genuinely_empty = not handle.read().strip()
            except (OSError, UnicodeError):
                genuinely_empty = False
        if allow_empty and genuinely_empty:
            contigs = ()
        else:
            raise IntegrityError(f"Cannot validate FASTA {fasta_path}: {error}") from error
    except (OSError, ValueError) as error:
        raise IntegrityError(f"Cannot validate FASTA {fasta_path}: {error}") from error
    fasta_ids = tuple(contig.sequence_id for contig in contigs)
    sequences = {contig.sequence_id: contig.sequence for contig in contigs}
    if len(sequences) != len(contigs):
        _fail(f"{fasta_path}: FASTA IDs are not unique")
    rows = _read_tsv(truth_path, TRUTH_COLUMNS)
    for row in rows:
        _validate_logical_row(
            row,
            SchemaName.TRUTH_ROW,
            integer_fields=("source_start", "source_end", "fragment_length"),
        )
    truth_ids = [row["sequence_id"] for row in rows]
    duplicates = sorted(identifier for identifier, count in Counter(truth_ids).items() if count > 1)
    if duplicates:
        _fail(f"{truth_path}: duplicate truth sequence IDs: {duplicates[:5]}")
    if tuple(truth_ids) != fasta_ids:
        missing_truth = sorted(set(sequences) - set(truth_ids))
        missing_fasta = sorted(set(truth_ids) - set(sequences))
        _fail(
            f"{split_dir.name}/{name}: FASTA and truth must have identical ordered IDs; "
            f"missing_truth={missing_truth[:5]}, missing_fasta={missing_fasta[:5]}"
        )
    required_view = expected_view or name
    semantic_partitions = allowed_partitions or frozenset({name})
    for line_number, row in enumerate(rows, start=2):
        identifier = row["sequence_id"]
        if re.fullmatch(r"frag-[0-9a-f]{32}", identifier) is None:
            _fail(f"{truth_path}:{line_number}: sequence_id is not an opaque CHIMERA ID")
        if row["view"] != required_view:
            _fail(f"{truth_path}:{line_number}: view must be {required_view!r}")
        if row["partition"] not in semantic_partitions:
            _fail(
                f"{truth_path}:{line_number}: semantic partition must be one of "
                f"{sorted(semantic_partitions)!r}"
            )
        try:
            Label(row["label"])
        except ValueError as error:
            raise IntegrityError(
                f"{truth_path}:{line_number}: label must be 'virus' or 'host'"
            ) from error
        start = _parse_int(
            row["source_start"], path=truth_path, line_number=line_number, field="source_start"
        )
        end = _parse_int(
            row["source_end"], path=truth_path, line_number=line_number, field="source_end"
        )
        length = _parse_int(
            row["fragment_length"],
            path=truth_path,
            line_number=line_number,
            field="fragment_length",
        )
        if start < 0 or end <= start or length <= 0 or end - start != length:
            _fail(f"{truth_path}:{line_number}: invalid zero-based half-open coordinates")
        if len(sequences[identifier]) != length:
            _fail(f"{truth_path}:{line_number}: fragment_length disagrees with FASTA sequence")
        if row["strand"] not in {"+", "-"}:
            _fail(f"{truth_path}:{line_number}: strand must be '+' or '-'")
        genome_id = row["source_genome_id"]
        source_id = row["source_sequence_id"]
        if not genome_id or not row["source_content_group_id"] or not source_id:
            _fail(f"{truth_path}:{line_number}: source identifiers must not be empty")
        source = sources.get(source_id)
        reference = references.get(genome_id)
        if source is None or reference is None or source.genome_id != genome_id:
            _fail(f"{truth_path}:{line_number}: source identifiers do not resolve")
        if row["label"] != source.label or row["label"] != reference.label:
            _fail(f"{truth_path}:{line_number}: label disagrees with source inventory")
        if row["source_accession_version"] != source.accession_version:
            _fail(f"{truth_path}:{line_number}: accession disagrees with source sequence")
        if row["release_date"] != source.release_date:
            _fail(f"{truth_path}:{line_number}: release_date disagrees with source sequence")
        if row["source_content_group_id"] != f"sha256:{reference.digest}":
            _fail(f"{truth_path}:{line_number}: source content group digest is incorrect")
        if length > len(source.sequence):
            _fail(f"{truth_path}:{line_number}: fragment exceeds its source sequence length")
        if source.topology == "linear":
            if end > len(source.sequence):
                _fail(f"{truth_path}:{line_number}: linear coordinates exceed source length")
            expected_coordinate_system = "0-based-half-open"
            expected_sequence = source.sequence[start:end]
        else:
            if start >= len(source.sequence) or end > start + len(source.sequence):
                _fail(f"{truth_path}:{line_number}: invalid circular source coordinates")
            expected_coordinate_system = "0-based-half-open-circular"
            expected_sequence = (source.sequence + source.sequence)[start:end]
        if row["coordinate_system"] != expected_coordinate_system:
            _fail(f"{truth_path}:{line_number}: coordinate_system disagrees with topology")
        if row["strand"] == "-":
            expected_sequence = reverse_complement(expected_sequence)
        if sequences[identifier] != expected_sequence:
            _fail(
                f"{truth_path}:{line_number}: FASTA sequence cannot be reconstructed "
                "from the declared source coordinates and strand"
            )
        if row["synthetic"] != "true":
            _fail(f"{truth_path}:{line_number}: synthetic must be 'true'")
    return _Partition(name=name, rows=rows, sequences=sequences, ordered_ids=fasta_ids)


def _compare_partition_stats(
    split_path: Path,
    split_manifest: Mapping[str, Any],
    partition: _Partition,
    *,
    check_truth_rows: bool = True,
) -> None:
    recorded = split_manifest.get(partition.name)
    if not isinstance(recorded, Mapping):
        _fail(f"{split_path}: missing {partition.name!r} statistics object")
    lengths = Counter(len(sequence) for sequence in partition.sequences.values())
    labels = Counter(row["label"] for row in partition.rows)
    label_lengths = Counter(
        f"{row['label']}:{len(partition.sequences[row['sequence_id']])}" for row in partition.rows
    )
    genomes = Counter(row["source_genome_id"] for row in partition.rows)
    total_bases = sum(length * count for length, count in lengths.items())
    gc_bases = sum(
        sequence.count("G") + sequence.count("C") for sequence in partition.sequences.values()
    )
    ambiguous_bases = sum(
        sum(base not in "ACGT" for base in sequence) for sequence in partition.sequences.values()
    )
    expected: dict[str, object] = {
        "records": len(partition.rows),
        "bases": total_bases,
        "gc_fraction": gc_bases / total_bases if total_bases else None,
        "ambiguous_fraction": ambiguous_bases / total_bases if total_bases else None,
        "records_by_label": dict(sorted(labels.items())),
        "records_by_label_and_length": dict(sorted(label_lengths.items())),
        "records_by_length": {str(length): count for length, count in sorted(lengths.items())},
        "records_by_genome": dict(sorted(genomes.items())),
        "source_genomes": len(genomes),
    }
    for key, value in expected.items():
        actual = recorded.get(key)
        if isinstance(value, float) and isinstance(actual, (int, float)):
            matches = math.isclose(float(actual), value, rel_tol=0.0, abs_tol=1e-15)
        else:
            matches = actual == value
        if not matches:
            _fail(
                f"{split_path}: recorded {partition.name}.{key} does not match "
                "independently recomputed data"
            )
    if check_truth_rows:
        recorded_truth_rows = split_manifest.get("truth_rows")
        if not isinstance(recorded_truth_rows, Mapping) or recorded_truth_rows.get(
            partition.name
        ) != len(partition.rows):
            _fail(f"{split_path}: recorded truth_rows.{partition.name} is incorrect")


def _validate_truth_assignments(
    kind: SplitKind,
    split_dir: Path,
    partitions: tuple[_Partition, _Partition],
    assignments: Mapping[str, Mapping[str, str]],
    references: Mapping[str, _Reference],
) -> None:
    train, test = partitions
    if train.ids & test.ids:
        _fail(f"{split_dir}: fragment IDs overlap train and test")
    for partition in partitions:
        labels = {row["label"] for row in partition.rows}
        if labels != {Label.VIRUS.value, Label.HOST.value}:
            _fail(f"{split_dir}/{partition.name}: both virus and host labels are required")
        for row in partition.rows:
            genome_id = row["source_genome_id"]
            assignment = assignments.get(genome_id)
            reference = references.get(genome_id)
            if assignment is None or reference is None:
                _fail(
                    f"{split_dir}/{partition.name}: truth names unassigned source group "
                    f"{genome_id!r}"
                )
            expected_partition = "both" if kind is SplitKind.RANDOM else partition.name
            if assignment["partition"] != expected_partition:
                _fail(
                    f"{split_dir}/{partition.name}: source group {genome_id!r} is assigned "
                    f"to {assignment['partition']!r}, expected {expected_partition!r}"
                )
            if row["label"] != assignment["label"] or row["label"] != reference.label:
                _fail(
                    f"{split_dir}/{partition.name}: truth label for {genome_id!r} "
                    "disagrees with its assignment/reference"
                )
            if row["source_sequence_id"] not in reference.sequence_ids:
                _fail(
                    f"{split_dir}/{partition.name}: truth names unknown source sequence "
                    f"{row['source_sequence_id']!r} for genome {genome_id!r}"
                )

    if kind is SplitKind.RANDOM:
        expected_sources = frozenset(references)
        if train.source_genomes != expected_sources or test.source_genomes != expected_sources:
            _fail(
                f"{split_dir}: every reference genome must occur in both random-fragment partitions"
            )
        train_strata = {(row["source_genome_id"], row["fragment_length"]) for row in train.rows}
        test_strata = {(row["source_genome_id"], row["fragment_length"]) for row in test.rows}
        if train_strata != test_strata:
            _fail(
                f"{split_dir}: every genome/length stratum must occur in both "
                "random-fragment partitions"
            )
    else:
        expected_train = frozenset(
            genome_id
            for genome_id, assignment in assignments.items()
            if assignment["partition"] == "train"
        )
        expected_test = frozenset(
            genome_id
            for genome_id, assignment in assignments.items()
            if assignment["partition"] == "test"
        )
        if train.source_genomes != expected_train:
            _fail(f"{split_dir}: training FASTA does not cover every assigned training genome")
        if test.source_genomes != expected_test:
            _fail(f"{split_dir}: test FASTA does not cover every assigned test genome")
        shared_genomes = train.source_genomes & test.source_genomes
        if shared_genomes:
            _fail(
                f"{split_dir}: source genomes leak across train/test: {sorted(shared_genomes)[:5]}"
            )
        train_hashes = {assignments[genome_id]["group_id"] for genome_id in train.source_genomes}
        test_hashes = {assignments[genome_id]["group_id"] for genome_id in test.source_genomes}
        if train_hashes & test_hashes:
            _fail(f"{split_dir}: canonical genome content leaks across train/test")
        if kind is SplitKind.SIMILARITY:
            train_fragment_hashes = {
                canonical_sequence_hash(sequence) for sequence in train.sequences.values()
            }
            test_fragment_hashes = {
                canonical_sequence_hash(sequence) for sequence in test.sequences.values()
            }
            if train_fragment_hashes & test_fragment_hashes:
                _fail(f"{split_dir}: exact fragment content leaks across strict train/test")

    content_labels: dict[str, set[str]] = {}
    for partition in partitions:
        for row in partition.rows:
            digest = canonical_sequence_hash(partition.sequences[row["sequence_id"]])
            content_labels.setdefault(digest, set()).add(row["label"])
    conflicts = [digest for digest, labels in content_labels.items() if len(labels) > 1]
    if conflicts:
        _fail(
            f"{split_dir}: {len(conflicts)} exact fragment sequence(s) carry contradictory "
            "virus/host labels"
        )


def _jsonable(value: object) -> object:
    if isinstance(value, date):
        return value.isoformat()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, SimilarityBands):
        return {"high": value.high, "moderate": value.moderate, "low": value.low}
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in sorted(value.items())}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return value


def _compare_assignment(
    path: Path,
    line_number: int,
    observed: Mapping[str, str],
    expected: SplitAssignment,
) -> None:
    expected_text = {
        "genome_id": expected.genome_id,
        "group_id": expected.group_id or expected.genome_id,
        "label": expected.label.value,
        "partition": expected.partition.value,
        "candidate_partition": (
            expected.candidate_partition.value if expected.candidate_partition else ""
        ),
        "reason": expected.reason,
        "release_date": expected.release_date.isoformat() if expected.release_date else "",
        "taxon": expected.taxon or "",
        "similarity_bin": expected.similarity_bin or "",
        "nearest_train_genome_id": expected.nearest_train_genome_id or "",
        "similarity_method": expected.similarity_method or "",
        "strict_gate_train_genome_id": expected.strict_gate_train_genome_id or "",
        "strict_gate_method": expected.strict_gate_method or "",
    }
    for field, value in expected_text.items():
        if observed[field] != value:
            _fail(
                f"{path}:{line_number}: {field} disagrees with independently "
                "recomputed split assignment"
            )
    expected_numbers = {
        "max_train_similarity": expected.max_train_similarity,
        "similarity_coverage": expected.similarity_coverage,
        "strict_gate_similarity": expected.strict_gate_similarity,
        "strict_gate_coverage": expected.strict_gate_coverage,
    }
    for numeric_field, numeric_value in expected_numbers.items():
        parsed = _parse_fraction(
            observed[numeric_field],
            path=path,
            line_number=line_number,
            field=numeric_field,
            allow_empty=True,
        )
        if parsed != numeric_value:
            _fail(
                f"{path}:{line_number}: {numeric_field} disagrees with independently "
                "recomputed similarity evidence"
            )


def _validate_recomputed_assignments(
    *,
    root: Path,
    split_dir: Path,
    kind: SplitKind,
    parameters: Mapping[str, Any],
    assignments: Mapping[str, Mapping[str, str]],
    genomes: tuple[Genome, ...],
    config: BenchmarkConfig,
) -> SplitPlan | None:
    """Re-run the canonical source-level split and compare every assignment field."""

    assignment_path = split_dir / "assignments.tsv"
    if kind is SplitKind.RANDOM:
        expected_parameters = {
            "diagnostic_only": True,
            "split_unit": "fragment",
            "test_fraction": config.test_fraction,
        }
        if dict(parameters) != expected_parameters:
            _fail(f"{split_dir / 'split.json'}: random parameters disagree with resolved config")
        for line_number, genome in enumerate(sorted(genomes, key=lambda item: item.genome_id), 2):
            observed = assignments[genome.genome_id]
            expected = SplitAssignment(
                genome_id=genome.genome_id,
                label=genome.label,
                partition=SplitPartition.TRAIN,
                reason="diagnostic_random_fragment_split",
                group_id=f"sha256:{genome.digest}",
                release_date=genome.metadata.release_date,
            )
            expected_text = {
                "group_id": expected.group_id or "",
                "label": expected.label.value,
                "partition": "both",
                "candidate_partition": "",
                "reason": expected.reason,
                "release_date": (
                    expected.release_date.isoformat() if expected.release_date else ""
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
            for field, value in expected_text.items():
                if observed[field] != value:
                    _fail(f"{assignment_path}:{line_number}: invalid random-fragment {field}")
        return None

    evidence_mode = parameters.get("similarity_evidence_mode")
    similarity_table: Path | None = None
    if kind is SplitKind.SIMILARITY and evidence_mode == "external-all-pairs":
        similarity_table = split_dir / "external-similarity.tsv"
        if not similarity_table.is_file():
            _fail(f"{split_dir}: external similarity evidence snapshot is missing")
    elif kind is SplitKind.SIMILARITY and evidence_mode != "built-in-recomputable":
        _fail(f"{split_dir / 'split.json'}: unknown similarity_evidence_mode")
    try:
        plan = build_split_plan(
            kind,
            genomes,
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
            similarity_table=similarity_table,
        )
    except ChimeraError as error:
        raise IntegrityError(f"{split_dir}: cannot reproduce split plan: {error}") from error

    expected_parameters = cast(dict[str, object], _jsonable(plan.parameters))
    if kind is SplitKind.SIMILARITY and similarity_table is not None:
        original = config.similarity_table.name if config.similarity_table is not None else None
        expected_parameters["similarity_source"] = original
        expected_parameters["similarity_table"] = original
    if dict(parameters) != expected_parameters:
        _fail(f"{split_dir / 'split.json'}: parameters disagree with resolved config/recomputation")
    expected_by_id = {assignment.genome_id: assignment for assignment in plan.assignments}
    for line_number, genome_id in enumerate(assignments, 2):
        _compare_assignment(
            assignment_path,
            line_number,
            assignments[genome_id],
            expected_by_id[genome_id],
        )
    return plan


def _validate_similarity_truth_row(
    split_dir: Path,
    row: Mapping[str, str],
    assignment: Mapping[str, str],
) -> None:
    genome_id = row["source_genome_id"]
    expected = {
        "partition": assignment["partition"],
        "similarity_bin": assignment["similarity_bin"],
        "max_train_similarity": assignment["max_train_similarity"],
        "nearest_train_genome_id": assignment["nearest_train_genome_id"],
    }
    for field, expected_value in expected.items():
        if row[field] != expected_value:
            _fail(
                f"{split_dir}: similarity truth field {field!r} for {genome_id!r} "
                "disagrees with assignments.tsv"
            )


def _validate_similarity_views(
    split_dir: Path,
    split_manifest: Mapping[str, Any],
    assignments: Mapping[str, Mapping[str, str]],
    strict_test: _Partition,
    sources: Mapping[str, _SourceSequence],
    references: Mapping[str, _Reference],
    expected_candidate: tuple[Fragment, ...],
) -> int:
    candidate = _read_partition(
        split_dir,
        "candidate_test",
        sources=sources,
        references=references,
        expected_view="candidate_test",
        allowed_partitions=frozenset({"test", "excluded"}),
    )
    _compare_expected_partition(split_dir, candidate, expected_candidate)
    expected_candidate_groups = {
        genome_id
        for genome_id, assignment in assignments.items()
        if assignment["candidate_partition"] == "test"
    }
    if candidate.source_genomes != expected_candidate_groups:
        _fail(f"{split_dir}: candidate_test does not contain exactly every proposed test genome")
    if {row["label"] for row in candidate.rows} != {
        Label.VIRUS.value,
        Label.HOST.value,
    }:
        _fail(f"{split_dir}/candidate_test: both virus and host labels are required")
    _compare_partition_stats(
        split_dir / "split.json",
        split_manifest,
        candidate,
        check_truth_rows=False,
    )
    for row in candidate.rows:
        _validate_similarity_truth_row(split_dir, row, assignments[row["source_genome_id"]])

    expected_strict_ids = {
        row["sequence_id"]
        for row in candidate.rows
        if assignments[row["source_genome_id"]]["partition"] == "test"
    }
    if strict_test.ids != expected_strict_ids:
        _fail(f"{split_dir}: strict test is not the correctly filtered candidate_test subset")
    if any(
        strict_test.sequences[identifier] != candidate.sequences[identifier]
        for identifier in strict_test.ids
    ):
        _fail(f"{split_dir}: strict test sequence differs from candidate_test")
    candidate_by_id = {row["sequence_id"]: row for row in candidate.rows}
    for row in strict_test.rows:
        _validate_similarity_truth_row(split_dir, row, assignments[row["source_genome_id"]])
        candidate_row = candidate_by_id[row["sequence_id"]]
        if any(row[field] != candidate_row[field] for field in TRUTH_COLUMNS if field != "view"):
            _fail(f"{split_dir}: strict test truth differs from candidate_test provenance")

    seen_ids: set[str] = set()
    strata_root = split_dir / "test_strata"
    auxiliary_records = len(candidate.rows)
    for similarity_bin in _SIMILARITY_BINS:
        stratum = _read_partition(
            strata_root,
            similarity_bin,
            sources=sources,
            references=references,
            expected_view=f"test_strata/{similarity_bin}",
            allowed_partitions=frozenset({"test", "excluded"}),
            allow_empty=True,
        )
        auxiliary_records += len(stratum.rows)
        overlap = seen_ids & set(stratum.ids)
        if overlap:
            _fail(f"{strata_root}: similarity strata overlap: {sorted(overlap)[:5]}")
        seen_ids.update(stratum.ids)
        expected_order = tuple(
            row["sequence_id"]
            for row in candidate.rows
            if assignments[row["source_genome_id"]]["similarity_bin"] == similarity_bin
        )
        if stratum.ordered_ids != expected_order:
            _fail(f"{strata_root}: {similarity_bin} order/content differs from candidate_test")
        for row in stratum.rows:
            assignment = assignments.get(row["source_genome_id"])
            if assignment is None or assignment["candidate_partition"] != "test":
                _fail(f"{strata_root}: stratum contains a non-candidate source genome")
            if assignment["similarity_bin"] != similarity_bin:
                _fail(
                    f"{strata_root}: genome {row['source_genome_id']!r} is in the wrong "
                    "similarity stratum"
                )
            _validate_similarity_truth_row(split_dir, row, assignment)
            stratum_candidate_row = candidate_by_id.get(row["sequence_id"])
            if stratum.sequences[row["sequence_id"]] != candidate.sequences.get(row["sequence_id"]):
                _fail(f"{strata_root}: stratum sequence differs from candidate_test")
            if stratum_candidate_row is None or any(
                row[field] != stratum_candidate_row[field]
                for field in TRUTH_COLUMNS
                if field != "view"
            ):
                _fail(f"{strata_root}: stratum truth differs from candidate_test provenance")
    if seen_ids != set(candidate.ids):
        missing = sorted(set(candidate.ids) - seen_ids)
        extra = sorted(seen_ids - set(candidate.ids))
        _fail(
            f"{strata_root}: strata must partition candidate_test exactly; "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )
    return auxiliary_records


def _validate_exclusions(
    path: Path,
    kind: SplitKind,
    assignments: Mapping[str, Mapping[str, str]],
    references: Mapping[str, _Reference],
) -> None:
    rows = _read_tsv(path, _EXCLUSION_COLUMNS)
    for row in rows:
        _validate_logical_row(row, SchemaName.EXCLUSION_ROW)
    by_id: dict[str, Mapping[str, str]] = {}
    for line_number, row in enumerate(rows, start=2):
        genome_id = row["genome_id"]
        if genome_id in by_id:
            _fail(f"{path}:{line_number}: duplicate excluded genome {genome_id!r}")
        assignment = assignments.get(genome_id)
        if assignment is None or assignment["partition"] != "excluded":
            _fail(f"{path}:{line_number}: exclusion is not backed by an excluded assignment")
        if row["split"] != kind.value or row["label"] != assignment["label"]:
            _fail(f"{path}:{line_number}: exclusion split/label disagrees with assignment")
        if row["reason"] != assignment["reason"]:
            _fail(f"{path}:{line_number}: exclusion reason disagrees with assignment")
        reference = references[genome_id]
        expected = {
            "duplicate_of": "",
            "source_sha256": reference.digest,
            "source_accession_version": reference.accession_version,
            "release_date": reference.release_date,
            "nearest_train_genome_id": assignment["nearest_train_genome_id"],
            "max_train_similarity": assignment["max_train_similarity"],
            "similarity_coverage": assignment["similarity_coverage"],
            "similarity_method": assignment["similarity_method"],
            "strict_gate_train_genome_id": assignment["strict_gate_train_genome_id"],
            "strict_gate_similarity": assignment["strict_gate_similarity"],
            "strict_gate_coverage": assignment["strict_gate_coverage"],
            "strict_gate_method": assignment["strict_gate_method"],
        }
        for field, expected_value in expected.items():
            if row[field] != expected_value:
                _fail(f"{path}:{line_number}: exclusion {field} disagrees with assignment/source")
        by_id[genome_id] = row
    expected_excluded = {
        genome_id
        for genome_id, assignment in assignments.items()
        if assignment["partition"] == "excluded"
    }
    if set(by_id) != expected_excluded:
        _fail(f"{path}: excluded.tsv must account for every excluded assignment exactly once")


def _semantic_fragment_order(
    fragments: tuple[Fragment, ...], *, seed: int, namespace: str
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


def _expected_fragment_views(
    kind: SplitKind,
    genomes: tuple[Genome, ...],
    plan: SplitPlan | None,
    config: BenchmarkConfig,
) -> tuple[tuple[Fragment, ...], tuple[Fragment, ...], tuple[Fragment, ...]]:
    by_id = {genome.genome_id: genome for genome in genomes}
    if kind is SplitKind.RANDOM:
        active = genomes
    else:
        assert plan is not None
        active_ids = set(plan.train_ids | plan.test_ids)
        if kind is SplitKind.SIMILARITY:
            active_ids.update(
                assignment.genome_id
                for assignment in plan.excluded
                if assignment.candidate_partition is SplitPartition.TEST
            )
        active = tuple(by_id[genome_id] for genome_id in sorted(active_ids))
    try:
        fragments = generate_fragments(
            active,
            fragment_lengths=config.fragment_lengths,
            fragments_per_genome=config.fragments_per_genome,
            seed=derive_seed(config.seed, "fragments", "shared-library-v1"),
            strand_mode=config.strand_mode,
            max_ambiguous_fraction=config.max_ambiguous_fraction,
        )
        if kind is SplitKind.RANDOM:
            train, test = split_fragments_random(
                fragments,
                test_fraction=config.test_fraction,
                seed=derive_seed(config.seed, "fragment-membership", kind.value),
            )
            candidate: tuple[Fragment, ...] = ()
        else:
            assert plan is not None
            train = tuple(
                fragment for fragment in fragments if fragment.genome_id in plan.train_ids
            )
            test = tuple(fragment for fragment in fragments if fragment.genome_id in plan.test_ids)
            candidate_ids = {
                assignment.genome_id
                for assignment in plan.assignments
                if assignment.candidate_partition is SplitPartition.TEST
            }
            candidate = tuple(
                fragment for fragment in fragments if fragment.genome_id in candidate_ids
            )
    except ChimeraError as error:
        raise IntegrityError(f"Cannot reproduce emitted fragments: {error}") from error
    return (
        _semantic_fragment_order(train, seed=config.seed, namespace=f"{kind.value}:train"),
        _semantic_fragment_order(test, seed=config.seed, namespace=f"{kind.value}:test"),
        (
            _semantic_fragment_order(
                candidate,
                seed=config.seed,
                namespace=f"{kind.value}:candidate-test",
            )
            if kind is SplitKind.SIMILARITY
            else ()
        ),
    )


def _compare_expected_partition(
    split_dir: Path,
    partition: _Partition,
    expected: tuple[Fragment, ...],
) -> None:
    expected_ids = tuple(fragment.fragment_id for fragment in expected)
    if partition.ordered_ids != expected_ids:
        _fail(
            f"{split_dir}/{partition.name}: emitted fragment IDs/order disagree with "
            "the resolved config and deterministic source sampler"
        )
    expected_by_id = {fragment.fragment_id: fragment for fragment in expected}
    for row in partition.rows:
        fragment = expected_by_id[row["sequence_id"]]
        expected_fields = {
            "label": fragment.label.value,
            "source_genome_id": fragment.genome_id,
            "source_sequence_id": fragment.sequence_id,
            "source_start": str(fragment.start),
            "source_end": str(fragment.end),
            "strand": fragment.strand,
            "fragment_length": str(fragment.length),
        }
        for field, value in expected_fields.items():
            if row[field] != value:
                _fail(
                    f"{split_dir}/{partition.name}: truth {field} for "
                    f"{fragment.fragment_id!r} disagrees with deterministic sampling"
                )


def _coordinate_segments(row: Mapping[str, str], source_length: int) -> tuple[tuple[int, int], ...]:
    start = int(row["source_start"])
    end = int(row["source_end"])
    if row["coordinate_system"] == "0-based-half-open" or end <= source_length:
        return ((start, end),)
    return ((start, source_length), (0, end - source_length))


def _coordinate_overlap_count(
    train: _Partition,
    test: _Partition,
    sources: Mapping[str, _SourceSequence],
) -> int:
    train_segments: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for row in train.rows:
        source = sources[row["source_sequence_id"]]
        key = (row["source_genome_id"], row["source_sequence_id"])
        train_segments.setdefault(key, []).extend(_coordinate_segments(row, len(source.sequence)))
    count = 0
    for row in test.rows:
        source = sources[row["source_sequence_id"]]
        key = (row["source_genome_id"], row["source_sequence_id"])
        test_segments = _coordinate_segments(row, len(source.sequence))
        if any(
            left_start < right_end and right_start < left_end
            for left_start, left_end in train_segments.get(key, ())
            for right_start, right_end in test_segments
        ):
            count += 1
    return count


def _validate_recorded_split_checks(
    split_path: Path,
    split_manifest: Mapping[str, Any],
    kind: SplitKind,
    train: _Partition,
    test: _Partition,
    assignments: Mapping[str, Mapping[str, str]],
    sources: Mapping[str, _SourceSequence],
) -> None:
    train_genomes = train.source_genomes
    test_genomes = test.source_genomes
    train_content = {canonical_sequence_hash(sequence) for sequence in train.sequences.values()}
    test_content = {canonical_sequence_hash(sequence) for sequence in test.sequences.values()}
    expected = {
        "status": "pass",
        "diagnostic_only": kind is SplitKind.RANDOM,
        "fragment_id_overlap": len(train.ids & test.ids),
        "exact_fragment_content_overlap": len(train_content & test_content),
        "test_fragments_with_coordinate_overlap": _coordinate_overlap_count(train, test, sources),
        "source_genome_overlap": len(train_genomes & test_genomes),
        "source_content_hash_overlap": len(
            {assignments[item]["group_id"] for item in train_genomes}
            & {assignments[item]["group_id"] for item in test_genomes}
        ),
        "train_source_genomes": len(train_genomes),
        "test_source_genomes": len(test_genomes),
    }
    if split_manifest.get("validation") != expected:
        _fail(f"{split_path}: validation summary disagrees with independent recomputation")


def _validate_split(
    root: Path,
    kind: SplitKind,
    embedded_manifest: Mapping[str, Any],
    references: Mapping[str, _Reference],
    sources: Mapping[str, _SourceSequence],
    genomes: tuple[Genome, ...],
    config: BenchmarkConfig,
) -> tuple[int, int, int, int]:
    split_dir = root / kind.directory_name
    required = (
        "assignments.tsv",
        "excluded.tsv",
        "split.json",
        "train.fasta.gz",
        "train.truth.tsv.gz",
        "test.fasta.gz",
        "test.truth.tsv.gz",
    )
    _require_files(split_dir, required)
    split_path = split_dir / "split.json"
    split_manifest = _read_json(split_path)
    validate_instance(split_manifest, SchemaName.SPLIT)
    if split_manifest != embedded_manifest:
        _fail(f"{split_path}: content disagrees with the copy embedded in manifest.json")
    if (
        split_manifest.get("schema") != _SPLIT_SCHEMA
        or split_manifest.get("protocol") != kind.value
        or split_manifest.get("protocol_id") != kind.value
    ):
        _fail(f"{split_path}: schema/protocol does not match directory {split_dir.name!r}")
    parameters = split_manifest.get("parameters")
    if not isinstance(parameters, Mapping):
        _fail(f"{split_path}: parameters must be an object")

    assignments = _read_assignments(split_dir / "assignments.tsv", kind, references)
    plan = _validate_recomputed_assignments(
        root=root,
        split_dir=split_dir,
        kind=kind,
        parameters=parameters,
        assignments=assignments,
        genomes=genomes,
        config=config,
    )
    _validate_exclusions(split_dir / "excluded.tsv", kind, assignments, references)
    train = _read_partition(
        split_dir,
        "train",
        sources=sources,
        references=references,
    )
    test = _read_partition(
        split_dir,
        "test",
        sources=sources,
        references=references,
    )
    expected_train, expected_test, expected_candidate = _expected_fragment_views(
        kind, genomes, plan, config
    )
    _compare_expected_partition(split_dir, train, expected_train)
    _compare_expected_partition(split_dir, test, expected_test)
    _validate_truth_assignments(kind, split_dir, (train, test), assignments, references)
    _compare_partition_stats(split_path, split_manifest, train)
    _compare_partition_stats(split_path, split_manifest, test)
    excluded = sum(row["partition"] == "excluded" for row in assignments.values())
    if split_manifest.get("excluded_genomes") != excluded:
        _fail(f"{split_path}: excluded_genomes count is incorrect")
    _validate_recorded_split_checks(
        split_path,
        split_manifest,
        kind,
        train,
        test,
        assignments,
        sources,
    )

    auxiliary_records = 0
    if kind is SplitKind.SIMILARITY:
        similarity_required = (
            "candidate_test.fasta.gz",
            "candidate_test.truth.tsv.gz",
            "test_strata/high_similarity.fasta.gz",
            "test_strata/high_similarity.truth.tsv.gz",
            "test_strata/moderate_similarity.fasta.gz",
            "test_strata/moderate_similarity.truth.tsv.gz",
            "test_strata/low_similarity.fasta.gz",
            "test_strata/low_similarity.truth.tsv.gz",
            "test_strata/distant_detectable.fasta.gz",
            "test_strata/distant_detectable.truth.tsv.gz",
            "test_strata/no_detectable_match.fasta.gz",
            "test_strata/no_detectable_match.truth.tsv.gz",
        )
        _require_files(split_dir, similarity_required)
        auxiliary_records = _validate_similarity_views(
            split_dir,
            split_manifest,
            assignments,
            test,
            sources,
            references,
            expected_candidate,
        )
    return len(train.rows), len(test.rows), len(assignments), auxiliary_records


def _read_resolved_config(root: Path, manifest: Mapping[str, Any]) -> BenchmarkConfig:
    path = root / "resolved-config.json"
    values = _read_json(path)
    validate_instance(values, SchemaName.RESOLVED_CONFIG)
    if "overwrite" in values:
        _fail(f"{path}: operational overwrite state must not be part of semantic provenance")
    try:
        config = BenchmarkConfig.from_mapping(values, base_dir=root)
    except (ChimeraError, TypeError, ValueError) as error:
        raise IntegrityError(f"{path}: invalid resolved configuration: {error}") from error
    split_manifests = manifest.get("splits")
    if not isinstance(split_manifests, Mapping):
        _fail(f"{root / 'manifest.json'}: splits must be an object")
    if set(split_manifests) != {kind.value for kind in config.splits}:
        _fail("resolved-config.json split set disagrees with manifest.json")
    randomness = manifest.get("randomness")
    if not isinstance(randomness, Mapping) or randomness.get("master_seed") != config.seed:
        _fail("resolved-config.json seed disagrees with manifest.json randomness")
    return config


def _validate_marker(root: Path) -> None:
    marker = root / ".chimera-bundle"
    try:
        value = marker.read_bytes()
    except OSError as error:
        raise IntegrityError(f"Cannot read bundle marker {marker}: {error}") from error
    expected = f"{_BUNDLE_SCHEMA}\n".encode("ascii")
    if value != expected:
        _fail(f"{marker}: marker must contain exactly {_BUNDLE_SCHEMA!r}")


def _validate_execution(root: Path, manifest: Mapping[str, Any]) -> None:
    path = root / "execution.json"
    execution = _read_json(path)
    required = {"started_at_utc", "finished_at_utc", "python", "platform", "command", "status"}
    if set(execution) != required or execution.get("status") != "complete":
        _fail(f"{path}: execution record is incomplete or has unexpected fields")
    timestamps: list[datetime] = []
    for field in ("started_at_utc", "finished_at_utc"):
        raw = execution.get(field)
        if not isinstance(raw, str):
            _fail(f"{path}: {field} must be an ISO timestamp")
        try:
            parsed = datetime.fromisoformat(raw)
        except ValueError as error:
            raise IntegrityError(f"{path}: {field} is not an ISO timestamp") from error
        if parsed.tzinfo is None or parsed.utcoffset() is None:
            _fail(f"{path}: {field} must include a UTC offset")
        timestamps.append(parsed)
    if timestamps[1] < timestamps[0]:
        _fail(f"{path}: finished_at_utc precedes started_at_utc")
    randomness = manifest.get("randomness")
    expected_python = randomness.get("python_version") if isinstance(randomness, Mapping) else None
    if execution.get("python") != expected_python:
        _fail(f"{path}: Python version disagrees with checked manifest provenance")
    command = execution.get("command")
    if not isinstance(command, list) or any(not isinstance(item, str) for item in command):
        _fail(f"{path}: command must be a JSON array of strings")
    if not isinstance(execution.get("platform"), str) or not execution["platform"]:
        _fail(f"{path}: platform must be a non-empty string")


def _validate_root_exclusions(
    root: Path,
    manifest: Mapping[str, Any],
    references: Mapping[str, _Reference],
) -> None:
    path = root / "excluded.tsv"
    rows = _read_tsv(path, _EXCLUSION_COLUMNS)
    for row in rows:
        _validate_logical_row(row, SchemaName.EXCLUSION_ROW)
    reference_manifest = manifest.get("references")
    if not isinstance(reference_manifest, Mapping) or reference_manifest.get(
        "preflight_exclusions"
    ) != len(rows):
        _fail("manifest references.preflight_exclusions disagrees with excluded.tsv")
    seen: set[str] = set()
    similarity_fields = {
        "nearest_train_genome_id",
        "max_train_similarity",
        "similarity_coverage",
        "similarity_method",
        "strict_gate_train_genome_id",
        "strict_gate_similarity",
        "strict_gate_coverage",
        "strict_gate_method",
    }
    for line_number, row in enumerate(rows, start=2):
        genome_id = row["genome_id"]
        if not genome_id or genome_id in seen or genome_id in references:
            _fail(f"{path}:{line_number}: excluded genome ID is empty, duplicate, or retained")
        seen.add(genome_id)
        try:
            label = Label(row["label"])
        except ValueError as error:
            raise IntegrityError(f"{path}:{line_number}: invalid exclusion label") from error
        duplicate_of = references.get(row["duplicate_of"])
        if (
            row["split"] != "reference_preflight"
            or row["reason"] != "same_class_content_duplicate"
            or duplicate_of is None
            or duplicate_of.label != label.value
        ):
            _fail(f"{path}:{line_number}: invalid preflight duplicate exclusion provenance")
        if row["source_sha256"] != duplicate_of.digest:
            _fail(f"{path}:{line_number}: duplicate digest disagrees with retained representative")
        if any(row[field] for field in similarity_fields):
            _fail(f"{path}:{line_number}: preflight exclusion carries split similarity evidence")
        if row["release_date"]:
            _parse_date(
                row["release_date"], path=path, line_number=line_number, field="release_date"
            )


def _validate_manifest_provenance(root: Path, manifest: Mapping[str, Any]) -> None:
    path = root / "manifest.json"
    tool = manifest.get("tool")
    if not isinstance(tool, Mapping) or tool.get("name") != "CHIMERA":
        _fail(f"{path}: invalid tool provenance")
    if tool.get("version") != __version__:
        _fail(f"{path}: bundle tool version is unsupported by this validator")
    if tool.get("software_content_sha256") != software_content_sha256():
        _fail(f"{path}: software content receipt disagrees with this CHIMERA installation")
    revision = tool.get("git_revision")
    if revision != "unknown" and (
        not isinstance(revision, str) or re.fullmatch(r"[0-9a-f]{40}", revision) is None
    ):
        _fail(f"{path}: invalid git_revision")
    git_dirty = tool.get("git_dirty")
    if git_dirty is not None and not isinstance(git_dirty, bool):
        _fail(f"{path}: git_dirty must be boolean or null")
    if revision == "unknown" and git_dirty is not None:
        _fail(f"{path}: git_dirty must be null when git_revision is unknown")
    randomness = manifest.get("randomness")
    if not isinstance(randomness, Mapping) or (
        randomness.get("algorithm")
        != "Python random.Random with semantic BLAKE2b-derived sub-seeds"
        or randomness.get("seed_derivation") != "chimera.seed.v1"
        or not isinstance(randomness.get("python_implementation"), str)
        or not isinstance(randomness.get("python_version"), str)
    ):
        _fail(f"{path}: invalid randomness provenance")
    references = manifest.get("references")
    inputs = references.get("inputs") if isinstance(references, Mapping) else None
    if not isinstance(inputs, list) or len(inputs) < 2:
        _fail(f"{path}: references.inputs must inventory source inputs")
    seen: set[tuple[str, str]] = set()
    for index, item in enumerate(inputs):
        if not isinstance(item, Mapping) or set(item) != {"content_id", "role", "sha256"}:
            _fail(f"{path}: references.inputs[{index}] is malformed")
        content_id, role, digest = (
            item.get("content_id"),
            item.get("role"),
            item.get("sha256"),
        )
        if (
            not isinstance(content_id, str)
            or not content_id
            or role not in {"reference_fasta", "reference_metadata", "external_similarity_table"}
            or not isinstance(digest, str)
            or re.fullmatch(r"[0-9a-f]{64}", digest) is None
            or content_id != f"sha256:{digest}"
            or (content_id, cast(str, role)) in seen
        ):
            _fail(f"{path}: references.inputs[{index}] has invalid provenance")
        seen.add((content_id, cast(str, role)))


def _validate_input_inventory_links(
    root: Path,
    manifest: Mapping[str, Any],
    config: BenchmarkConfig,
    sources: Mapping[str, _SourceSequence],
) -> None:
    references = manifest["references"]
    assert isinstance(references, Mapping)
    inputs = references["inputs"]
    assert isinstance(inputs, list)
    by_role: dict[str, set[str]] = {}
    for item in inputs:
        assert isinstance(item, Mapping)
        by_role.setdefault(cast(str, item["role"]), set()).add(cast(str, item["content_id"]))
    virus_sources = {
        source.source_input_id for source in sources.values() if source.label == Label.VIRUS.value
    }
    host_sources = {
        source.source_input_id for source in sources.values() if source.label == Label.HOST.value
    }
    configured_virus = {path.name for path in config.virus_paths}
    configured_host = {path.name for path in config.host_paths}
    if not virus_sources <= configured_virus:
        _fail(f"{root / 'resolved-config.json'}: viral input IDs disagree with sequences.tsv")
    if not host_sources <= configured_host:
        _fail(f"{root / 'resolved-config.json'}: host input IDs disagree with sequences.tsv")
    if configured_virus & configured_host:
        _fail(f"{root / 'resolved-config.json'}: an input is assigned to both biological labels")
    if by_role.get("reference_fasta", set()) != configured_virus | configured_host:
        _fail(f"{root / 'manifest.json'}: reference input receipts disagree with config")
    expected_metadata = {config.metadata_path.name} if config.metadata_path is not None else set()
    if by_role.get("reference_metadata", set()) != expected_metadata:
        _fail(f"{root / 'manifest.json'}: metadata input provenance disagrees with config")
    expected_similarity = (
        {config.similarity_table.name} if config.similarity_table is not None else set()
    )
    if by_role.get("external_similarity_table", set()) != expected_similarity:
        _fail(f"{root / 'manifest.json'}: similarity input provenance disagrees with config")


def _validate_embedded_schemas(root: Path) -> None:
    """Require a byte-semantic copy of every schema needed to read the bundle."""

    schema_dir = root / "schemas"
    if not schema_dir.is_dir():
        _fail(f"{schema_dir}: embedded schema directory is missing")
    expected = {schema_filename(name) for name in JSON_SCHEMA_NAMES}
    actual = {path.name for path in schema_dir.iterdir()}
    if actual != expected:
        _fail(
            f"{schema_dir}: embedded schema inventory is incomplete or contains extras; "
            f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}"
        )
    for schema_name in JSON_SCHEMA_NAMES:
        path = schema_dir / schema_filename(schema_name)
        embedded = _read_json(path)
        if embedded != load_schema(schema_name):
            _fail(f"{path}: embedded schema differs from CHIMERA {__version__}")


def _validate_exact_directory(path: Path, expected: set[str] | frozenset[str]) -> None:
    """Reject undeclared files and directories from a versioned bundle location."""

    try:
        actual = {child.name for child in path.iterdir()}
    except OSError as error:
        raise IntegrityError(f"Cannot enumerate bundle directory {path}: {error}") from error
    if actual != expected:
        _fail(
            f"{path}: non-canonical bundle layout; "
            f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}"
        )


def _validate_bundle_layout(root: Path, config: BenchmarkConfig) -> None:
    """Require the exact v1 tree, including conditional similarity evidence."""

    split_directories = {kind.directory_name for kind in config.splits}
    root_entries = set(_ROOT_FILES) | {"schemas"} | split_directories
    _validate_exact_directory(root, root_entries)
    _validate_exact_directory(
        root / "schemas",
        {schema_filename(name) for name in JSON_SCHEMA_NAMES},
    )

    for kind in config.splits:
        split_dir = root / kind.directory_name
        expected = set(_PRIMARY_SPLIT_FILES)
        if kind is SplitKind.SIMILARITY:
            expected.update(
                {
                    "candidate_test.fasta.gz",
                    "candidate_test.truth.tsv.gz",
                    "test_strata",
                }
            )
            if config.similarity_table is not None:
                expected.add("external-similarity.tsv")
        _validate_exact_directory(split_dir, expected)
        if kind is SplitKind.SIMILARITY:
            strata_files = {
                f"{similarity_bin}.{suffix}"
                for similarity_bin in _SIMILARITY_BINS
                for suffix in ("fasta.gz", "truth.tsv.gz")
            }
            _validate_exact_directory(split_dir / "test_strata", strata_files)


def validate_bundle(root: Path | str) -> ValidationReport:
    """Independently validate a published CHIMERA benchmark directory.

    Checksums are verified first.  Every split named by the root manifest must
    then have the canonical on-disk layout, exact FASTA/truth ID agreement,
    consistent assignments and labels, valid coordinates, both classes in both
    partitions, recomputed counts, and its protocol-specific leakage invariant.

    Raises:
        IntegrityError: If any checksum, file, schema, or scientific invariant
            fails.  Recorded manifest pass flags are intentionally ignored.
    """

    requested_root = Path(root).expanduser()
    if requested_root.is_symlink():
        _fail(f"Bundle root must not be a symbolic link: {requested_root}")
    try:
        root_path = requested_root.resolve(strict=True)
    except OSError as error:
        raise IntegrityError(f"Cannot access bundle root {requested_root}: {error}") from error
    if not root_path.is_dir():
        _fail(f"Bundle root is not a directory: {root_path}")
    _require_files(root_path, _ROOT_FILES)
    checksum_result = verify_checksums(root_path)
    checksums_verified = checksum_result.get("files_checked")
    if isinstance(checksums_verified, bool) or not isinstance(checksums_verified, int):
        _fail("Checksum verifier returned an invalid file count")

    _validate_embedded_schemas(root_path)
    manifest_path = root_path / "manifest.json"
    manifest = _read_json(manifest_path)
    validate_instance(manifest, SchemaName.BUNDLE)
    if manifest.get("schema") != _BUNDLE_SCHEMA:
        _fail(f"{manifest_path}: unsupported bundle schema {manifest.get('schema')!r}")
    split_manifests = manifest.get("splits")
    if not isinstance(split_manifests, Mapping) or not split_manifests:
        _fail(f"{manifest_path}: splits must be a non-empty object")
    _validate_marker(root_path)
    _validate_manifest_provenance(root_path, manifest)
    config = _read_resolved_config(root_path, manifest)
    _validate_bundle_layout(root_path, config)
    _validate_execution(root_path, manifest)
    references = _read_references(root_path, manifest)
    sources = _read_source_sequences(root_path, references)
    genomes = _reconstruct_genomes(references, sources)
    _validate_input_inventory_links(root_path, manifest, config, sources)
    _validate_root_exclusions(root_path, manifest, references)
    observed_split_directories = {
        child.name
        for child in root_path.iterdir()
        if child.is_dir() and child.name in {kind.directory_name for kind in SplitKind}
    }
    expected_split_directories = {kind.directory_name for kind in config.splits}
    if observed_split_directories != expected_split_directories:
        _fail("On-disk protocol directories disagree with the resolved split set")

    split_counts: list[tuple[str, int, int]] = []
    primary_fasta_records = 0
    primary_truth_rows = 0
    auxiliary_fasta_records = 0
    auxiliary_truth_rows = 0
    assignment_rows = 0
    checks = [
        "checksums",
        "bundle_manifest",
        "resolved_config",
        "execution_record",
        "reference_inventory",
        "source_sequence_authentication",
        "preflight_exclusions",
    ]
    for raw_kind, embedded in split_manifests.items():
        try:
            kind = SplitKind(str(raw_kind))
        except ValueError as error:
            raise IntegrityError(f"{manifest_path}: unknown split protocol {raw_kind!r}") from error
        if not isinstance(embedded, Mapping):
            _fail(f"{manifest_path}: split {kind.value!r} must be an object")
        train_count, test_count, assignment_count, auxiliary_count = _validate_split(
            root_path,
            kind,
            embedded,
            references,
            sources,
            genomes,
            config,
        )
        split_counts.append((kind.value, train_count, test_count))
        primary_fasta_records += train_count + test_count
        primary_truth_rows += train_count + test_count
        auxiliary_fasta_records += auxiliary_count
        auxiliary_truth_rows += auxiliary_count
        assignment_rows += assignment_count
        checks.extend(
            (
                f"{kind.value}:fasta_truth",
                f"{kind.value}:assignments",
                f"{kind.value}:counts",
                f"{kind.value}:protocol_invariants",
            )
        )

    split_counts.sort(key=lambda item: list(SplitKind).index(SplitKind(item[0])))
    return ValidationReport(
        root=root_path,
        checksums_verified=checksums_verified,
        split_counts=tuple(split_counts),
        primary_fasta_records_verified=primary_fasta_records,
        primary_truth_rows_verified=primary_truth_rows,
        auxiliary_fasta_records_verified=auxiliary_fasta_records,
        auxiliary_truth_rows_verified=auxiliary_truth_rows,
        assignment_rows_verified=assignment_rows,
        checks=tuple(checks),
    )


__all__ = ["ValidationReport", "validate_bundle"]
