"""Reference discovery, metadata joining, grouping, and contamination checks."""

from __future__ import annotations

import csv
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from datetime import date
from pathlib import Path

from .config import DuplicatePolicy
from .errors import InputError
from .fasta import read_fasta
from .models import (
    Contig,
    Genome,
    GenomeMetadata,
    Label,
    Topology,
    deterministic_topology_agnostic_genome_hash,
)

TAXONOMY_RANKS = (
    "realm",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)

METADATA_REQUIRED_COLUMNS = ("sequence_id",)
METADATA_RECOMMENDED_COLUMNS = (
    "genome_id",
    "label",
    "accession_version",
    "release_date",
    "topology",
    *TAXONOMY_RANKS,
)


@dataclass(frozen=True, slots=True)
class MetadataRecord:
    """Normalized metadata for one FASTA record/segment."""

    sequence_id: str
    genome_id: str
    label: Label | None
    accession_version: str | None
    release_date: date | None
    topology: Topology
    taxonomy: tuple[tuple[str, str], ...]
    extra: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True)
class ReferenceExclusion:
    """A reference removed before split planning, always with an audit reason."""

    genome_id: str
    label: Label
    reason: str
    duplicate_of: str | None = None
    source_sha256: str = ""
    accession_version: str | None = None
    release_date: date | None = None


@dataclass(frozen=True, slots=True)
class ReferenceCatalog:
    """Validated, content-deduplicated source collection."""

    genomes: tuple[Genome, ...]
    source_files: tuple[Path, ...]
    exclusions: tuple[ReferenceExclusion, ...] = ()
    warnings: tuple[str, ...] = ()

    def by_id(self) -> dict[str, Genome]:
        """Return genomes indexed by globally unique stable group ID."""

        return {genome.genome_id: genome for genome in self.genomes}


def _normalized_header(fieldnames: Iterable[str | None]) -> dict[str, str]:
    result: dict[str, str] = {}
    for raw in fieldnames:
        if raw is None:
            continue
        canonical = raw.strip().lower().replace("-", "_").replace(" ", "_")
        if canonical in result:
            raise InputError(f"Metadata contains duplicate column {canonical!r}")
        result[canonical] = raw
    return result


def _parse_release_date(raw: str, *, path: Path, line_number: int) -> date | None:
    value = raw.strip()
    if not value:
        return None
    try:
        return date.fromisoformat(value)
    except ValueError as exc:
        raise InputError(
            f"{path}:{line_number}: release_date must use ISO YYYY-MM-DD, got {value!r}"
        ) from exc


def read_metadata(path: Path) -> dict[str, MetadataRecord]:
    """Read sequence-level TSV/CSV metadata with strict schema validation.

    ``release_date`` means the accession's first public release date. Ambiguous
    legacy date names are rejected so deposition/creation timestamps cannot be
    silently relabeled as public availability.
    """

    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    try:
        handle = path.open(newline="", encoding="utf-8-sig")
    except OSError as exc:
        raise InputError(f"Cannot read metadata table {path}: {exc}") from exc
    records: dict[str, MetadataRecord] = {}
    with handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        header = _normalized_header(reader.fieldnames or ())
        if "sequence_id" not in header:
            raise InputError(f"Metadata {path} is missing required column 'sequence_id'")
        date_columns = [
            name for name in ("release_date", "deposited_at", "create_date") if name in header
        ]
        legacy_date_columns = [name for name in ("deposited_at", "create_date") if name in header]
        if legacy_date_columns:
            raise InputError(
                f"Metadata {path} uses ambiguous legacy date column(s) {legacy_date_columns}; "
                "provide a verified first-public date in release_date"
            )
        if len(date_columns) > 1:
            raise InputError(
                f"Metadata {path} has ambiguous date columns {date_columns}; use only release_date"
            )
        for line_number, raw_row in enumerate(reader, start=2):
            row = {
                canonical: (raw_row.get(original) or "").strip()
                for canonical, original in header.items()
            }
            sequence_id = row["sequence_id"]
            if not sequence_id:
                raise InputError(f"{path}:{line_number}: sequence_id must not be empty")
            if sequence_id in records:
                raise InputError(
                    f"{path}:{line_number}: duplicate metadata sequence_id {sequence_id!r}"
                )
            genome_id = row.get("genome_id") or sequence_id
            raw_label = row.get("label", "")
            try:
                label = Label(raw_label.lower()) if raw_label else None
            except ValueError as exc:
                raise InputError(
                    f"{path}:{line_number}: label must be 'virus' or 'host', got {raw_label!r}"
                ) from exc
            date_column = date_columns[0] if date_columns else ""
            release_date = _parse_release_date(
                row.get(date_column, ""), path=path, line_number=line_number
            )
            topology_value = (row.get("topology") or "linear").lower()
            if topology_value == "linear":
                topology: Topology = "linear"
            elif topology_value == "circular":
                topology = "circular"
            else:
                raise InputError(f"{path}:{line_number}: topology must be 'linear' or 'circular'")
            taxonomy = tuple(
                (rank, row.get(rank) or row.get(f"tax_{rank}", ""))
                for rank in TAXONOMY_RANKS
                if row.get(rank) or row.get(f"tax_{rank}")
            )
            consumed = {
                "sequence_id",
                "genome_id",
                "label",
                "accession_version",
                "release_date",
                "deposited_at",
                "create_date",
                "topology",
                *TAXONOMY_RANKS,
                *(f"tax_{rank}" for rank in TAXONOMY_RANKS),
            }
            extra = tuple(sorted((key, value) for key, value in row.items() if key not in consumed))
            try:
                records[sequence_id] = MetadataRecord(
                    sequence_id=sequence_id,
                    genome_id=genome_id,
                    label=label,
                    accession_version=row.get("accession_version") or None,
                    release_date=release_date,
                    topology=topology,
                    taxonomy=taxonomy,
                    extra=extra,
                )
            except (TypeError, ValueError) as exc:
                raise InputError(f"{path}:{line_number}: {exc}") from exc
    if not records:
        raise InputError(f"Metadata table {path} has no data rows")
    return records


def _metadata_for_group(
    genome_id: str,
    records: list[MetadataRecord],
) -> GenomeMetadata:
    dates = sorted(record.release_date for record in records if record.release_date is not None)
    # A segmented/grouped genome becomes publicly usable only when every
    # segment has a known release date; its effective date is then the latest.
    # A partial date set is unknown, never an invitation to impute from the
    # available segments.
    all_segment_dates_known = len(dates) == len(records)
    group_release_date = dates[-1] if dates and all_segment_dates_known else None
    taxonomy_by_rank: dict[str, set[str]] = defaultdict(set)
    for record in records:
        for rank, value in record.taxonomy:
            taxonomy_by_rank[rank].add(value)
    conflicts = {rank: values for rank, values in taxonomy_by_rank.items() if len(values) > 1}
    if conflicts:
        detail = "; ".join(
            f"{rank}={sorted(values)!r}" for rank, values in sorted(conflicts.items())
        )
        raise InputError(f"Genome group {genome_id!r} has inconsistent taxonomy: {detail}")
    taxonomy = tuple(
        (rank, next(iter(values))) for rank, values in sorted(taxonomy_by_rank.items())
    )
    accession_values = {
        record.accession_version for record in records if record.accession_version is not None
    }
    accession_version = next(iter(accession_values)) if len(accession_values) == 1 else None
    extra_common: dict[str, str] = {}
    if records:
        segment_extras = [dict(record.extra) for record in records]
        common_keys = set(segment_extras[0])
        for extra in segment_extras[1:]:
            common_keys.intersection_update(extra)
        extra_common = {
            key: segment_extras[0][key]
            for key in common_keys
            if all(extra[key] == segment_extras[0][key] for extra in segment_extras[1:])
        }
    if len(records) > 1:
        extra_common["grouped_sequence_count"] = str(len(records))
        if all_segment_dates_known:
            extra_common["group_release_date_policy"] = "latest_segment_release"
        else:
            extra_common["group_release_date_policy"] = "unknown_incomplete_segment_dates"
            extra_common["missing_release_date_sequences"] = ",".join(
                sorted(record.sequence_id for record in records if record.release_date is None)
            )
    return GenomeMetadata(
        release_date=group_release_date,
        taxonomy=taxonomy,
        accession_version=accession_version,
        extra=tuple(sorted(extra_common.items())),
    )


def _load_labeled_contigs(
    paths: tuple[Path, ...],
    label: Label,
    metadata: Mapping[str, MetadataRecord] | None,
) -> tuple[list[tuple[Contig, MetadataRecord]], tuple[Path, ...]]:
    contigs = read_fasta(paths)
    if metadata is None:
        records_by_file: dict[Path, int] = defaultdict(int)
        for contig in contigs:
            if contig.source_path is None:
                raise InputError(f"FASTA sequence {contig.sequence_id!r} lacks a source path")
            records_by_file[contig.source_path] += 1
        multi_record_files = sorted(
            (path for path, count in records_by_file.items() if count > 1),
            key=lambda path: path.as_posix(),
        )
        if multi_record_files:
            raise InputError(
                "Metadata is required for multi-record FASTA because CHIMERA cannot infer "
                "whether records are contigs/segments of one genome or independent genomes: "
                + ", ".join(str(path) for path in multi_record_files[:5])
            )
    discovered_files = tuple(
        sorted({contig.source_path for contig in contigs if contig.source_path})
    )
    joined: list[tuple[Contig, MetadataRecord]] = []
    for contig in contigs:
        record = metadata.get(contig.sequence_id) if metadata is not None else None
        if metadata is not None and record is None:
            raise InputError(
                f"FASTA sequence {contig.sequence_id!r} has no row in the metadata table"
            )
        if record is None:
            record = MetadataRecord(
                sequence_id=contig.sequence_id,
                genome_id=contig.sequence_id,
                label=label,
                accession_version=contig.sequence_id,
                release_date=None,
                topology="linear",
                taxonomy=(),
                extra=(),
            )
        if record.label is not None and record.label != label:
            raise InputError(
                f"Metadata labels sequence {contig.sequence_id!r} as {record.label.value}, "
                f"but it was supplied via --{label.value}"
            )
        if record.label is None:
            record = MetadataRecord(
                sequence_id=record.sequence_id,
                genome_id=record.genome_id,
                label=label,
                accession_version=record.accession_version,
                release_date=record.release_date,
                topology=record.topology,
                taxonomy=record.taxonomy,
                extra=record.extra,
            )
        joined.append(
            (
                replace(
                    contig,
                    accession_version=record.accession_version,
                    release_date=record.release_date,
                    topology=record.topology,
                    taxonomy=record.taxonomy,
                    metadata_extra=record.extra,
                ),
                record,
            )
        )
    return joined, discovered_files


def _merge_same_class_duplicate_metadata(duplicates: list[Genome]) -> Genome:
    """Choose an earliest-evidence representative and preserve merge provenance."""

    taxonomy_signatures = {genome.metadata.taxonomy for genome in duplicates}
    if len(taxonomy_signatures) != 1:
        identifiers = ", ".join(sorted(genome.genome_id for genome in duplicates))
        raise InputError(
            "Same-class identical content has conflicting taxonomy and cannot be dropped: "
            f"{identifiers}"
        )
    representative = min(
        duplicates,
        key=lambda genome: (
            genome.metadata.release_date is None,
            genome.metadata.release_date or date.max,
            genome.genome_id,
        ),
    )
    known_dates = sorted(
        genome.metadata.release_date
        for genome in duplicates
        if genome.metadata.release_date is not None
    )
    merged_extra = dict(representative.metadata.extra)
    merged_extra.update(
        {
            "chimera_duplicate_genome_ids": ",".join(
                sorted(genome.genome_id for genome in duplicates)
            ),
            "chimera_duplicate_accession_versions": ",".join(
                sorted(
                    {
                        contig.accession_version
                        for genome in duplicates
                        for contig in genome.contigs
                        if contig.accession_version is not None
                    }
                )
            ),
            "chimera_duplicate_release_dates": ",".join(
                f"{genome.genome_id}:{genome.metadata.release_date.isoformat() if genome.metadata.release_date else 'unknown'}"
                for genome in sorted(duplicates, key=lambda item: item.genome_id)
            ),
            "chimera_duplicate_source_sequence_digests": ",".join(
                sorted({contig.digest for genome in duplicates for contig in genome.contigs})
            ),
            "chimera_duplicate_merge_policy": "earliest_known_public_release",
        }
    )
    return Genome(
        genome_id=representative.genome_id,
        label=representative.label,
        contigs=representative.contigs,
        metadata=GenomeMetadata(
            release_date=known_dates[0] if known_dates else None,
            taxonomy=representative.metadata.taxonomy,
            accession_version=representative.metadata.accession_version,
            extra=tuple(sorted(merged_extra.items())),
        ),
    )


def load_reference_catalog(
    virus_paths: tuple[Path, ...],
    host_paths: tuple[Path, ...],
    *,
    metadata_path: Path | None = None,
    duplicate_policy: DuplicatePolicy = DuplicatePolicy.ERROR,
) -> ReferenceCatalog:
    """Load, group, and content-audit virus/host references.

    Identical or reverse-complement-equivalent content with contradictory
    labels is always fatal.  Same-class duplicate genomes either fail or are
    deterministically removed according to ``duplicate_policy``.
    """

    metadata = read_metadata(metadata_path) if metadata_path is not None else None
    virus_rows, virus_files = _load_labeled_contigs(virus_paths, Label.VIRUS, metadata)
    host_rows, host_files = _load_labeled_contigs(host_paths, Label.HOST, metadata)
    joined = virus_rows + host_rows
    sequence_id_counts: dict[str, int] = defaultdict(int)
    for contig, _ in joined:
        sequence_id_counts[contig.sequence_id] += 1
    duplicate_sequence_ids = sorted(
        sequence_id for sequence_id, count in sequence_id_counts.items() if count > 1
    )
    if duplicate_sequence_ids:
        raise InputError(
            "Sequence IDs must be globally unique across virus and host inputs; duplicates: "
            + ", ".join(duplicate_sequence_ids[:5])
        )
    observed_ids = {contig.sequence_id for contig, _ in joined}
    if metadata is not None:
        unused = sorted(set(metadata) - observed_ids)
        if unused:
            preview = ", ".join(unused[:5])
            suffix = " …" if len(unused) > 5 else ""
            raise InputError(
                f"Metadata contains {len(unused)} sequence(s) absent from FASTA input: "
                f"{preview}{suffix}"
            )

    grouped_contigs: dict[str, list[Contig]] = defaultdict(list)
    grouped_records: dict[str, list[MetadataRecord]] = defaultdict(list)
    grouped_labels: dict[str, set[Label]] = defaultdict(set)
    for contig, record in joined:
        grouped_contigs[record.genome_id].append(contig)
        grouped_records[record.genome_id].append(record)
        assert record.label is not None
        grouped_labels[record.genome_id].add(record.label)
    conflicting_groups = {
        genome_id: labels for genome_id, labels in grouped_labels.items() if len(labels) > 1
    }
    if conflicting_groups:
        genome_id = sorted(conflicting_groups)[0]
        labels = sorted(label.value for label in conflicting_groups[genome_id])
        raise InputError(f"Genome group ID {genome_id!r} occurs in contradictory classes: {labels}")

    genomes: list[Genome] = []
    for genome_id in sorted(grouped_contigs):
        label = next(iter(grouped_labels[genome_id]))
        try:
            genomes.append(
                Genome(
                    genome_id=genome_id,
                    label=label,
                    contigs=tuple(
                        sorted(grouped_contigs[genome_id], key=lambda item: item.sequence_id)
                    ),
                    metadata=_metadata_for_group(genome_id, grouped_records[genome_id]),
                )
            )
        except (TypeError, ValueError) as exc:
            raise InputError(f"Invalid genome group {genome_id!r}: {exc}") from exc

    raw_digest_groups: dict[str, list[Genome]] = defaultdict(list)
    digest_groups: dict[str, list[Genome]] = defaultdict(list)
    for genome in genomes:
        raw_digest = deterministic_topology_agnostic_genome_hash(genome.contigs)
        raw_digest_groups[raw_digest].append(genome)
        digest_groups[genome.digest].append(genome)
    for raw_digest, exact_matches in sorted(raw_digest_groups.items()):
        exact_labels = {genome.label for genome in exact_matches}
        if len(exact_labels) > 1:
            ids = ", ".join(
                f"{genome.genome_id} ({genome.label.value}, "
                f"{','.join(sorted({contig.topology for contig in genome.contigs}))})"
                for genome in sorted(exact_matches, key=lambda item: item.genome_id)
            )
            raise InputError(
                "Cross-class content conflict: topology-agnostic identical/"
                "reverse-complement-equivalent genome content is labeled both virus and host: "
                f"{ids}; audit digest {raw_digest}"
            )
    kept: list[Genome] = []
    exclusions: list[ReferenceExclusion] = []
    for digest, duplicates in sorted(digest_groups.items()):
        duplicate_labels = {genome.label for genome in duplicates}
        if len(duplicate_labels) > 1:
            ids = ", ".join(
                f"{genome.genome_id} ({genome.label.value})"
                for genome in sorted(duplicates, key=lambda item: item.genome_id)
            )
            raise InputError(
                "Cross-class content conflict: identical/reverse-complement-equivalent "
                f"genome digest {digest} is labeled both virus and host: {ids}"
            )
        ordered = sorted(duplicates, key=lambda item: item.genome_id)
        if len(ordered) > 1 and duplicate_policy == DuplicatePolicy.ERROR:
            ids = ", ".join(genome.genome_id for genome in ordered)
            raise InputError(
                f"Same-class duplicate genome content detected ({digest}): {ids}; "
                "review the references or use --duplicate-policy drop"
            )
        representative = (
            _merge_same_class_duplicate_metadata(ordered) if len(ordered) > 1 else ordered[0]
        )
        kept.append(representative)
        exclusions.extend(
            ReferenceExclusion(
                genome_id=duplicate.genome_id,
                label=duplicate.label,
                reason="same_class_content_duplicate",
                duplicate_of=representative.genome_id,
                source_sha256=duplicate.digest,
                accession_version=duplicate.metadata.accession_version,
                release_date=duplicate.metadata.release_date,
            )
            for duplicate in ordered
            if duplicate.genome_id != representative.genome_id
        )
    kept.sort(key=lambda genome: (genome.label.value, genome.genome_id))
    source_files = tuple(sorted(set(virus_files + host_files), key=lambda path: path.as_posix()))
    return ReferenceCatalog(
        genomes=tuple(kept),
        source_files=source_files,
        exclusions=tuple(sorted(exclusions, key=lambda item: item.genome_id)),
    )
