"""Leakage-aware, deterministic genome split planning.

The functions in this module assign *source genomes*, never generated
fragments.  Generation therefore happens only after a plan has been written
and audited.  Exact sequence duplicates (including reverse-complement-only
representations) share a content group and cannot cross train/test boundaries.

Temporal plans use the first public ``release_date``.  Unless an immutable
historical database snapshot is declared, such a plan is correctly described
as a release-date-filtered retrospective split rather than a prospective one.
"""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from datetime import date
from enum import StrEnum
from pathlib import Path
from typing import TypeAlias, cast

from .config import MissingMetadataPolicy, SimilarityBands, SplitKind
from .errors import ConfigurationError, InputError, IntegrityError
from .models import Genome, Label
from .similarity import SimilarityHit, best_train_matches, read_similarity_table

ParameterValue: TypeAlias = (
    str | int | float | bool | date | Path | SimilarityBands | tuple[object, ...] | None
)


@dataclass(frozen=True, slots=True)
class FrozenParameters(Mapping[str, ParameterValue]):
    """Small immutable mapping used for resolved split provenance."""

    _entries: tuple[tuple[str, ParameterValue], ...] = ()

    def __post_init__(self) -> None:
        raw_entries: Iterable[tuple[str, object]]
        raw_entries = self._entries.items() if isinstance(self._entries, Mapping) else self._entries

        normalized: list[tuple[str, ParameterValue]] = []
        seen: set[str] = set()
        for raw_key, raw_value in raw_entries:
            key = str(raw_key)
            if not key or key in seen:
                raise ValueError(f"split parameter key {key!r} is empty or duplicated")
            seen.add(key)
            value: object = raw_value
            if isinstance(value, list):
                value = tuple(value)
            elif isinstance(value, set):
                value = tuple(sorted(value, key=str))
            if not isinstance(
                value,
                (str, int, float, bool, date, Path, SimilarityBands, tuple, type(None)),
            ):
                raise TypeError(
                    f"split parameter {key!r} has unsupported value type {type(value).__name__}"
                )
            normalized.append((key, cast(ParameterValue, value)))
        object.__setattr__(self, "_entries", tuple(sorted(normalized)))

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> FrozenParameters:
        """Copy *values* into an immutable, key-sorted mapping."""

        entries = tuple((key, cast(ParameterValue, value)) for key, value in values.items())
        return cls(entries)

    def __getitem__(self, key: str) -> ParameterValue:
        for candidate, value in self._entries:
            if candidate == key:
                return value
        raise KeyError(key)

    def __iter__(self) -> Iterator[str]:
        return (key for key, _ in self._entries)

    def __len__(self) -> int:
        return len(self._entries)


def _coerce_parameters(value: object) -> FrozenParameters:
    if isinstance(value, FrozenParameters):
        return value
    if isinstance(value, Mapping):
        return FrozenParameters.from_mapping(cast(Mapping[str, object], value))
    raise TypeError("parameters must be a mapping")


class SplitPartition(StrEnum):
    """A source genome's explicit disposition in a split plan."""

    TRAIN = "train"
    TEST = "test"
    EXCLUDED = "excluded"


_SIMILARITY_BINS = frozenset(
    {
        "high_similarity",
        "moderate_similarity",
        "low_similarity",
        "distant_detectable",
        "no_detectable_match",
    }
)
_MIN_GROUPS_PER_LABEL = 2
_ERROR_PREVIEW_COUNT = 5
_DEFAULT_SIMILARITY_BANDS = SimilarityBands()


@dataclass(frozen=True, slots=True)
class SplitAssignment:
    """One immutable, auditable source-genome assignment.

    ``candidate_partition`` records the genome-disjoint proposal used by the
    similarity protocol.  It remains ``test`` when the strict similarity gate
    changes the final partition to ``excluded``.
    """

    genome_id: str
    label: Label
    partition: SplitPartition
    reason: str
    group_id: str | None = None
    candidate_partition: SplitPartition | None = None
    release_date: date | None = None
    taxon: str | None = None
    similarity_bin: str | None = None
    nearest_train_genome_id: str | None = None
    max_train_similarity: float | None = None
    similarity_coverage: float | None = None
    similarity_method: str | None = None
    strict_gate_train_genome_id: str | None = None
    strict_gate_similarity: float | None = None
    strict_gate_coverage: float | None = None
    strict_gate_method: str | None = None

    def __post_init__(self) -> None:
        _normalize_assignment_enums(self)
        _normalize_assignment_text(self)
        _validate_assignment_metadata(self)


def _normalize_assignment_enums(assignment: SplitAssignment) -> None:
    try:
        object.__setattr__(assignment, "label", Label(assignment.label))
    except (TypeError, ValueError) as exc:
        raise ValueError("assignment label must be 'virus' or 'host'") from exc
    try:
        object.__setattr__(
            assignment,
            "partition",
            SplitPartition(assignment.partition),
        )
    except (TypeError, ValueError) as exc:
        raise ValueError("partition must be 'train', 'test', or 'excluded'") from exc
    if assignment.candidate_partition is None:
        return
    try:
        candidate = SplitPartition(assignment.candidate_partition)
    except (TypeError, ValueError) as exc:
        raise ValueError("candidate_partition must be 'train' or 'test'") from exc
    if candidate is SplitPartition.EXCLUDED:
        raise ValueError("candidate_partition cannot be 'excluded'")
    object.__setattr__(assignment, "candidate_partition", candidate)


def _normalize_assignment_text(assignment: SplitAssignment) -> None:
    if not isinstance(assignment.genome_id, str) or not assignment.genome_id:
        raise ValueError("genome_id must be a non-empty string")
    if not isinstance(assignment.reason, str) or not assignment.reason.strip():
        raise ValueError("assignment reason must be a non-empty string")
    if "\n" in assignment.reason or "\r" in assignment.reason or "\x00" in assignment.reason:
        raise ValueError("assignment reason must be one line")
    object.__setattr__(assignment, "reason", assignment.reason.strip())
    if assignment.group_id is None:
        object.__setattr__(assignment, "group_id", assignment.genome_id)
    elif not isinstance(assignment.group_id, str) or not assignment.group_id:
        raise ValueError("group_id must be a non-empty string")
    if assignment.taxon is not None:
        normalized_taxon = assignment.taxon.strip()
        if not normalized_taxon:
            raise ValueError("taxon cannot be blank")
        object.__setattr__(assignment, "taxon", normalized_taxon)


def _validate_assignment_metadata(assignment: SplitAssignment) -> None:
    if assignment.release_date is not None and type(assignment.release_date) is not date:
        raise TypeError("release_date must be datetime.date or None")
    if assignment.similarity_bin is not None and assignment.similarity_bin not in _SIMILARITY_BINS:
        raise ValueError(f"unknown similarity bin {assignment.similarity_bin!r}")
    for field_name in (
        "max_train_similarity",
        "similarity_coverage",
        "strict_gate_similarity",
        "strict_gate_coverage",
    ):
        value = getattr(assignment, field_name)
        if value is not None and (not math.isfinite(value) or not 0.0 <= value <= 1.0):
            raise ValueError(f"{field_name} must be finite and within [0, 1]")
    if assignment.nearest_train_genome_id is not None and not assignment.nearest_train_genome_id:
        raise ValueError("nearest_train_genome_id cannot be blank")
    if assignment.similarity_method is not None and not assignment.similarity_method.strip():
        raise ValueError("similarity_method cannot be blank")
    if assignment.strict_gate_similarity is None:
        if any(
            value is not None
            for value in (
                assignment.strict_gate_train_genome_id,
                assignment.strict_gate_coverage,
                assignment.strict_gate_method,
            )
        ):
            raise ValueError("strict gate provenance requires strict_gate_similarity")
    elif not assignment.strict_gate_train_genome_id or not assignment.strict_gate_method:
        raise ValueError("strict_gate_similarity requires a training genome ID and method")


@dataclass(frozen=True, slots=True)
class SplitPlan:
    """A complete immutable source partition and its resolved parameters."""

    kind: SplitKind
    assignments: tuple[SplitAssignment, ...]
    seed: int
    parameters: FrozenParameters = FrozenParameters()

    def __post_init__(self) -> None:
        try:
            object.__setattr__(self, "kind", SplitKind(self.kind))
        except (TypeError, ValueError) as exc:
            raise ValueError("kind must be a supported SplitKind") from exc
        if isinstance(self.seed, bool) or not isinstance(self.seed, int):
            raise TypeError("seed must be an integer")
        if isinstance(cast(object, self.assignments), (str, bytes)):
            raise TypeError("assignments must be an iterable of SplitAssignment objects")
        assignments = tuple(self.assignments)
        if not assignments or not all(
            isinstance(assignment, SplitAssignment) for assignment in assignments
        ):
            raise TypeError("assignments must contain SplitAssignment objects")
        object.__setattr__(
            self,
            "assignments",
            tuple(sorted(assignments, key=lambda assignment: assignment.genome_id)),
        )
        object.__setattr__(self, "parameters", _coerce_parameters(cast(object, self.parameters)))
        self.validate()

    @property
    def train(self) -> tuple[SplitAssignment, ...]:
        """Training assignments in stable genome-ID order."""

        return tuple(a for a in self.assignments if a.partition is SplitPartition.TRAIN)

    @property
    def test(self) -> tuple[SplitAssignment, ...]:
        """Retained test assignments in stable genome-ID order."""

        return tuple(a for a in self.assignments if a.partition is SplitPartition.TEST)

    @property
    def excluded(self) -> tuple[SplitAssignment, ...]:
        """Explicitly excluded assignments in stable genome-ID order."""

        return tuple(a for a in self.assignments if a.partition is SplitPartition.EXCLUDED)

    @property
    def train_ids(self) -> frozenset[str]:
        return frozenset(assignment.genome_id for assignment in self.train)

    @property
    def test_ids(self) -> frozenset[str]:
        return frozenset(assignment.genome_id for assignment in self.test)

    @property
    def excluded_ids(self) -> frozenset[str]:
        return frozenset(assignment.genome_id for assignment in self.excluded)

    def assignment_for(self, genome_id: str) -> SplitAssignment:
        """Return one assignment or raise ``KeyError`` for an unknown genome."""

        for assignment in self.assignments:
            if assignment.genome_id == genome_id:
                return assignment
        raise KeyError(genome_id)

    def validate(self) -> None:
        """Run all partition and protocol-specific leakage validators."""

        _validate_partition(self)
        if self.kind is SplitKind.TEMPORAL:
            _validate_temporal(self)
        elif self.kind is SplitKind.TAXONOMY:
            _validate_taxonomy(self)
        elif self.kind is SplitKind.SIMILARITY:
            _validate_similarity(self)


@dataclass(frozen=True, slots=True)
class _GenomeGroup:
    group_id: str
    label: Label
    genomes: tuple[Genome, ...]


def _semantic_rank(seed: int, namespace: str, value: str) -> bytes:
    """Return a process- and ordering-independent rank for a semantic key."""

    payload = f"v1\0{seed}\0{namespace}\0{value}".encode()
    return hashlib.blake2b(payload, digest_size=16, person=b"CHIMERA-split-v1").digest()


def _validate_fraction(test_fraction: float) -> None:
    if not math.isfinite(test_fraction) or not 0.0 < test_fraction < 1.0:
        raise ConfigurationError("test_fraction must be finite and strictly between 0 and 1")


def _coerce_missing_policy(value: MissingMetadataPolicy | str) -> MissingMetadataPolicy:
    try:
        return MissingMetadataPolicy(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError("missing_metadata must be 'error' or 'exclude'") from exc


def _prepare_groups(genomes: Sequence[Genome]) -> tuple[_GenomeGroup, ...]:
    if isinstance(cast(object, genomes), (str, bytes)):
        raise TypeError("genomes must be a sequence of Genome objects")
    source = tuple(genomes)
    if not source:
        raise ConfigurationError("At least one source genome is required")
    if not all(isinstance(genome, Genome) for genome in source):
        raise TypeError("genomes must contain only Genome objects")
    genome_ids = [genome.genome_id for genome in source]
    duplicated_ids = sorted(
        genome_id for genome_id in set(genome_ids) if genome_ids.count(genome_id) > 1
    )
    if duplicated_ids:
        raise InputError("Duplicate genome_id values: " + ", ".join(duplicated_ids))

    by_digest: dict[str, list[Genome]] = defaultdict(list)
    for genome in source:
        by_digest[genome.digest].append(genome)
    groups: list[_GenomeGroup] = []
    for digest, members in sorted(by_digest.items()):
        labels = {member.label for member in members}
        if len(labels) != 1:
            ids = ", ".join(sorted(member.genome_id for member in members))
            raise InputError(
                "Identical sequence content has conflicting virus/host labels; "
                f"quarantine or correct these genomes: {ids}"
            )
        groups.append(
            _GenomeGroup(
                group_id=f"sha256:{digest}",
                label=next(iter(labels)),
                genomes=tuple(sorted(members, key=lambda genome: genome.genome_id)),
            )
        )
    return tuple(groups)


def _group_lookup(groups: Sequence[_GenomeGroup]) -> dict[str, _GenomeGroup]:
    return {genome.genome_id: group for group in groups for genome in group.genomes}


def _partition_one_label(
    groups: Sequence[_GenomeGroup],
    *,
    test_fraction: float,
    seed: int,
    namespace: str,
    label: Label,
) -> dict[str, SplitPartition]:
    selected = [group for group in groups if group.label is label]
    if len(selected) < _MIN_GROUPS_PER_LABEL:
        raise ConfigurationError(
            f"{label.value} requires at least two independent genome content groups "
            "to place at least one genome in both train and test"
        )
    ranked = sorted(
        selected,
        key=lambda group: (
            _semantic_rank(seed, namespace, group.group_id),
            group.group_id,
        ),
    )
    target_genomes = sum(len(group.genomes) for group in ranked) * test_fraction
    cumulative = 0
    candidates: list[tuple[float, bytes, int]] = []
    for index, group in enumerate(ranked[:-1], start=1):
        cumulative += len(group.genomes)
        candidates.append(
            (
                abs(cumulative - target_genomes),
                _semantic_rank(seed, f"{namespace}:cut", str(index)),
                index,
            )
        )
    test_group_count = min(candidates)[2]
    test_groups = {group.group_id for group in ranked[:test_group_count]}
    return {
        genome.genome_id: (
            SplitPartition.TEST if group.group_id in test_groups else SplitPartition.TRAIN
        )
        for group in selected
        for genome in group.genomes
    }


def _stratified_partitions(
    groups: Sequence[_GenomeGroup],
    *,
    test_fraction: float,
    seed: int,
    namespace: str,
) -> dict[str, SplitPartition]:
    _validate_fraction(test_fraction)
    result: dict[str, SplitPartition] = {}
    for label in Label:
        result.update(
            _partition_one_label(
                groups,
                test_fraction=test_fraction,
                seed=seed,
                namespace=f"{namespace}:{label.value}",
                label=label,
            )
        )
    return result


def _release_date(genome: Genome) -> date | None:
    """Read canonical release_date, with a transition fallback for old models."""

    value = getattr(genome.metadata, "release_date", None)
    if value is None and not hasattr(genome.metadata, "release_date"):
        value = getattr(genome.metadata, "deposited_at", None)
    return value


def _taxonomy_value(genome: Genome, rank: str) -> str | None:
    return genome.metadata.taxon(rank)


def _parameters(values: Mapping[str, object]) -> FrozenParameters:
    return FrozenParameters.from_mapping(values)


def genome_holdout(
    genomes: Sequence[Genome],
    *,
    test_fraction: float = 0.20,
    seed: int = 42,
) -> SplitPlan:
    """Create Test 2B: label-stratified, content-group-disjoint holdout."""

    groups = _prepare_groups(genomes)
    partitions = _stratified_partitions(
        groups,
        test_fraction=test_fraction,
        seed=seed,
        namespace="genome-holdout",
    )
    lookup = _group_lookup(groups)
    assignments = tuple(
        SplitAssignment(
            genome_id=genome.genome_id,
            label=genome.label,
            partition=partitions[genome.genome_id],
            reason="label_stratified_genome_holdout",
            group_id=lookup[genome.genome_id].group_id,
            release_date=_release_date(genome),
        )
        for genome in sorted(genomes, key=lambda item: item.genome_id)
    )
    return SplitPlan(
        kind=SplitKind.GENOME,
        assignments=assignments,
        seed=seed,
        parameters=_parameters(
            {
                "grouping": "canonical_topology_aware_genome_sha256_v2",
                "test_fraction": test_fraction,
            }
        ),
    )


def _choose_temporal_cutoff(
    groups: Sequence[_GenomeGroup],
    *,
    dates: Mapping[str, date],
    test_fraction: float,
    seed: int,
) -> date:
    effective_dates = {
        group.group_id: min(
            dates[genome.genome_id] for genome in group.genomes if genome.genome_id in dates
        )
        for group in groups
        if any(genome.genome_id in dates for genome in group.genomes)
    }
    possible = sorted(set(effective_dates.values()))
    viable: list[tuple[float, float, bytes, date]] = []
    for cutoff in possible:
        fractions: list[float] = []
        all_test = 0
        all_count = 0
        is_viable = True
        for label in Label:
            label_groups = [
                group
                for group in groups
                if group.label is label and group.group_id in effective_dates
            ]
            train_count = sum(
                len([genome for genome in group.genomes if genome.genome_id in dates])
                for group in label_groups
                if effective_dates[group.group_id] <= cutoff
            )
            test_count = sum(
                len([genome for genome in group.genomes if genome.genome_id in dates])
                for group in label_groups
                if effective_dates[group.group_id] > cutoff
            )
            if train_count == 0 or test_count == 0:
                is_viable = False
                break
            fractions.append(test_count / (train_count + test_count))
            all_test += test_count
            all_count += train_count + test_count
        if is_viable:
            viable.append(
                (
                    sum(abs(fraction - test_fraction) for fraction in fractions),
                    abs((all_test / all_count) - test_fraction),
                    _semantic_rank(seed, "temporal-cutoff", cutoff.isoformat()),
                    cutoff,
                )
            )
    if not viable:
        raise ConfigurationError(
            "No release-date cutoff can place at least one virus and one host genome "
            "on both sides; provide more date-diverse data or an explicit viable cutoff"
        )
    return min(viable)[3]


def temporal_holdout(
    genomes: Sequence[Genome],
    *,
    test_fraction: float = 0.20,
    seed: int = 42,
    release_date: date | None = None,
    temporal_cutoff: date | None = None,
    missing_metadata: MissingMetadataPolicy | str = MissingMetadataPolicy.ERROR,
    historical_snapshot: str | None = None,
) -> SplitPlan:
    """Create Test 2D with an inclusive training release-date cutoff."""

    _validate_fraction(test_fraction)
    policy = _coerce_missing_policy(missing_metadata)
    if release_date is not None and temporal_cutoff is not None and release_date != temporal_cutoff:
        raise ConfigurationError("release_date and temporal_cutoff disagree")
    requested_cutoff = release_date if release_date is not None else temporal_cutoff
    if requested_cutoff is not None and type(requested_cutoff) is not date:
        raise ConfigurationError("release_date must be datetime.date or None")
    if historical_snapshot is not None and not historical_snapshot.strip():
        raise ConfigurationError("historical_snapshot cannot be blank")

    groups = _prepare_groups(genomes)
    lookup = _group_lookup(groups)
    dates = {
        genome.genome_id: value
        for genome in genomes
        if (value := _release_date(genome)) is not None
    }
    missing_ids = sorted(genome.genome_id for genome in genomes if genome.genome_id not in dates)
    if missing_ids and policy is MissingMetadataPolicy.ERROR:
        preview = ", ".join(missing_ids[:_ERROR_PREVIEW_COUNT])
        suffix = " …" if len(missing_ids) > _ERROR_PREVIEW_COUNT else ""
        raise ConfigurationError(
            f"{len(missing_ids)} genome(s) lack release_date: {preview}{suffix}; "
            "supply metadata or use missing_metadata='exclude'"
        )
    if not dates:
        raise ConfigurationError("No genomes have release_date metadata")

    cutoff = requested_cutoff or _choose_temporal_cutoff(
        groups,
        dates=dates,
        test_fraction=test_fraction,
        seed=seed,
    )
    effective_dates = {
        group.group_id: min(
            dates[genome.genome_id] for genome in group.genomes if genome.genome_id in dates
        )
        for group in groups
        if any(genome.genome_id in dates for genome in group.genomes)
    }

    assignments: list[SplitAssignment] = []
    for genome in sorted(genomes, key=lambda item: item.genome_id):
        value = dates.get(genome.genome_id)
        group = lookup[genome.genome_id]
        if value is None:
            partition = SplitPartition.EXCLUDED
            reason = "missing_release_date"
        elif effective_dates[group.group_id] <= cutoff:
            partition = SplitPartition.TRAIN
            reason = (
                "release_date_on_or_before_cutoff"
                if value <= cutoff
                else "sequence_group_available_on_or_before_cutoff"
            )
        else:
            partition = SplitPartition.TEST
            reason = "release_date_after_cutoff"
        assignments.append(
            SplitAssignment(
                genome_id=genome.genome_id,
                label=genome.label,
                partition=partition,
                reason=reason,
                group_id=group.group_id,
                release_date=value,
            )
        )

    semantics = (
        "historical-snapshot prospective"
        if historical_snapshot is not None
        else "release-date-filtered retrospective"
    )
    return SplitPlan(
        kind=SplitKind.TEMPORAL,
        assignments=tuple(assignments),
        seed=seed,
        parameters=_parameters(
            {
                "cutoff_selection": "explicit" if requested_cutoff is not None else "automatic",
                "grouping": "canonical_topology_aware_genome_sha256_v2",
                "historical_snapshot": historical_snapshot,
                "missing_metadata": policy.value,
                "release_date": cutoff,
                "temporal_semantics": semantics,
                "test_fraction": test_fraction,
            }
        ),
    )


def _normalize_requested_taxa(values: Iterable[str]) -> tuple[str, ...]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        normalized = value.strip()
        if not normalized:
            raise ConfigurationError("holdout_taxa cannot contain blank names")
        folded = normalized.casefold()
        if folded not in seen:
            seen.add(folded)
            result.append(normalized)
    return tuple(result)


def _collect_viral_taxa(
    genomes: Sequence[Genome],
    *,
    rank: str,
    policy: MissingMetadataPolicy,
) -> dict[str, str]:
    taxa = {
        genome.genome_id: value
        for genome in genomes
        if genome.label is Label.VIRUS and (value := _taxonomy_value(genome, rank)) is not None
    }
    missing = sorted(
        genome.genome_id
        for genome in genomes
        if genome.label is Label.VIRUS and genome.genome_id not in taxa
    )
    if missing and policy is MissingMetadataPolicy.ERROR:
        preview = ", ".join(missing[:_ERROR_PREVIEW_COUNT])
        suffix = " …" if len(missing) > _ERROR_PREVIEW_COUNT else ""
        raise ConfigurationError(
            f"{len(missing)} viral genome(s) lack taxonomy rank {rank!r}: "
            f"{preview}{suffix}; supply metadata or use missing_metadata='exclude'"
        )
    return taxa


def _validate_group_taxonomy(
    groups: Sequence[_GenomeGroup],
    *,
    taxa: Mapping[str, str],
    rank: str,
) -> None:
    for group in groups:
        if group.label is not Label.VIRUS:
            continue
        group_taxa = {
            taxa[genome.genome_id].casefold()
            for genome in group.genomes
            if genome.genome_id in taxa
        }
        if len(group_taxa) > 1:
            ids = ", ".join(genome.genome_id for genome in group.genomes)
            raise InputError(
                f"Identical viral sequence content has conflicting {rank} taxonomy: {ids}"
            )


def _select_holdout_taxa(
    taxa: Mapping[str, str],
    *,
    requested: tuple[str, ...],
    rank: str,
    count: int,
    seed: int,
) -> tuple[tuple[str, ...], str]:
    available: dict[str, str] = {}
    for taxon in sorted(set(taxa.values()), key=lambda value: (value.casefold(), value)):
        available.setdefault(taxon.casefold(), taxon)
    if len(available) < _MIN_GROUPS_PER_LABEL:
        raise ConfigurationError(
            f"Taxonomic holdout at rank {rank!r} requires at least two represented viral taxa"
        )
    if requested:
        unknown = [taxon for taxon in requested if taxon.casefold() not in available]
        if unknown:
            raise ConfigurationError(
                f"Requested {rank} holdout taxa are absent: {', '.join(unknown)}"
            )
        selected = tuple(available[taxon.casefold()] for taxon in requested)
        if len(selected) >= len(available):
            raise ConfigurationError(
                "Explicit holdout_taxa must leave at least one represented viral taxon in training"
            )
        return selected, "explicit"
    if count >= len(available):
        raise ConfigurationError(
            "auto_holdout_count must leave at least one represented viral taxon in training"
        )
    ranked = sorted(
        available.values(),
        key=lambda taxon: (
            _semantic_rank(seed, f"taxonomy:{rank}", taxon.casefold()),
            taxon.casefold(),
        ),
    )
    return tuple(ranked[:count]), "automatic"


def _make_taxonomic_assignments(
    genomes: Sequence[Genome],
    *,
    lookup: Mapping[str, _GenomeGroup],
    taxa: Mapping[str, str],
    selected: tuple[str, ...],
    rank: str,
    host_partitions: Mapping[str, SplitPartition],
) -> tuple[SplitAssignment, ...]:
    selected_folded = {taxon.casefold() for taxon in selected}
    assignments: list[SplitAssignment] = []
    for genome in sorted(genomes, key=lambda item: item.genome_id):
        taxon = taxa.get(genome.genome_id) if genome.label is Label.VIRUS else None
        if genome.label is Label.HOST:
            partition = host_partitions[genome.genome_id]
            reason = "host_label_stratified_genome_holdout"
        elif taxon is None:
            partition = SplitPartition.EXCLUDED
            reason = f"missing_taxonomy_{rank}"
        elif taxon.casefold() in selected_folded:
            partition = SplitPartition.TEST
            reason = "selected_taxon_holdout"
        else:
            partition = SplitPartition.TRAIN
            reason = "non_holdout_taxon"
        assignments.append(
            SplitAssignment(
                genome_id=genome.genome_id,
                label=genome.label,
                partition=partition,
                reason=reason,
                group_id=lookup[genome.genome_id].group_id,
                release_date=_release_date(genome),
                taxon=taxon,
            )
        )
    return tuple(assignments)


def taxonomic_holdout(
    genomes: Sequence[Genome],
    *,
    test_fraction: float = 0.20,
    seed: int = 42,
    taxonomy_rank: str = "family",
    holdout_taxa: Iterable[str] = (),
    auto_holdout_count: int = 1,
    missing_metadata: MissingMetadataPolicy | str = MissingMetadataPolicy.ERROR,
) -> SplitPlan:
    """Create Test 2E by holding complete viral taxa out of training."""

    _validate_fraction(test_fraction)
    policy = _coerce_missing_policy(missing_metadata)
    rank = taxonomy_rank.strip().lower()
    if not rank or any(character.isspace() for character in rank):
        raise ConfigurationError("taxonomy_rank must be one non-empty rank name")
    if isinstance(auto_holdout_count, bool) or auto_holdout_count < 1:
        raise ConfigurationError("auto_holdout_count must be at least 1")
    requested = _normalize_requested_taxa(holdout_taxa)

    groups = _prepare_groups(genomes)
    lookup = _group_lookup(groups)
    taxa = _collect_viral_taxa(genomes, rank=rank, policy=policy)
    _validate_group_taxonomy(groups, taxa=taxa, rank=rank)
    selected, selection = _select_holdout_taxa(
        taxa,
        requested=requested,
        rank=rank,
        count=auto_holdout_count,
        seed=seed,
    )

    host_partitions = _partition_one_label(
        groups,
        test_fraction=test_fraction,
        seed=seed,
        namespace="taxonomy:host",
        label=Label.HOST,
    )
    assignments = _make_taxonomic_assignments(
        genomes,
        lookup=lookup,
        taxa=taxa,
        selected=selected,
        rank=rank,
        host_partitions=host_partitions,
    )

    return SplitPlan(
        kind=SplitKind.TAXONOMY,
        assignments=assignments,
        seed=seed,
        parameters=_parameters(
            {
                "grouping": "canonical_topology_aware_genome_sha256_v2",
                "holdout_selection": selection,
                "holdout_taxa": selected,
                "missing_metadata": policy.value,
                "taxonomy_rank": rank,
                "test_fraction": test_fraction,
            }
        ),
    )


def _validate_similarity_parameters(
    *,
    similarity_k: int,
    sketch_size: int,
    max_train_similarity: float,
    min_similarity_coverage: float,
    similarity_bands: SimilarityBands,
) -> None:
    if isinstance(similarity_k, bool) or similarity_k < 1:
        raise ConfigurationError("similarity_k must be a positive integer")
    if isinstance(sketch_size, bool) or sketch_size < 1:
        raise ConfigurationError("sketch_size must be a positive integer")
    if not math.isfinite(max_train_similarity) or not 0.0 <= max_train_similarity <= 1.0:
        raise ConfigurationError("max_train_similarity must be finite and within [0, 1]")
    if not math.isfinite(min_similarity_coverage) or not 0.0 <= min_similarity_coverage <= 1.0:
        raise ConfigurationError("min_similarity_coverage must be finite and within [0, 1]")
    if not isinstance(similarity_bands, SimilarityBands):
        raise TypeError("similarity_bands must be a SimilarityBands object")
    if max_train_similarity < similarity_bands.high:
        raise ConfigurationError(
            "max_train_similarity must be at least the high novelty-band threshold"
        )


def _similarity_candidate_assignment(
    genome: Genome,
    *,
    group_id: str,
    hit: SimilarityHit,
    max_train_similarity: float,
    min_similarity_coverage: float,
    similarity_bands: SimilarityBands,
) -> SplitAssignment:
    if hit.query_genome_id != genome.genome_id:
        raise IntegrityError(
            f"Similarity hit for {genome.genome_id!r} names query {hit.query_genome_id!r}"
        )
    gate_reference = hit.strict_gate_reference_genome_id
    gate_similarity = hit.strict_gate_similarity
    gate_coverage = hit.strict_gate_coverage
    gate_method = hit.strict_gate_method
    maximum_above_threshold = hit.similarity is not None and hit.similarity > max_train_similarity
    maximum_has_coverage = hit.coverage is None or hit.coverage >= min_similarity_coverage
    if gate_similarity is None and maximum_above_threshold and maximum_has_coverage:
        gate_reference = hit.reference_genome_id
        gate_similarity = hit.similarity
        gate_coverage = hit.coverage
        gate_method = hit.method
    exclude = gate_similarity is not None
    if exclude:
        reason = "similarity_above_strict_identity_and_coverage_gate"
    elif maximum_above_threshold:
        reason = "similarity_above_identity_below_minimum_coverage"
    else:
        reason = "similarity_within_strict_identity_threshold"
    return SplitAssignment(
        genome_id=genome.genome_id,
        label=genome.label,
        partition=SplitPartition.EXCLUDED if exclude else SplitPartition.TEST,
        candidate_partition=SplitPartition.TEST,
        reason=reason,
        group_id=group_id,
        release_date=_release_date(genome),
        similarity_bin=similarity_bands.classify(hit.similarity),
        nearest_train_genome_id=hit.reference_genome_id,
        max_train_similarity=hit.similarity,
        similarity_coverage=hit.coverage,
        similarity_method=hit.method,
        strict_gate_train_genome_id=gate_reference,
        strict_gate_similarity=gate_similarity,
        strict_gate_coverage=gate_coverage,
        strict_gate_method=gate_method,
    )


def similarity_filtered_holdout(
    genomes: Sequence[Genome],
    *,
    test_fraction: float = 0.20,
    seed: int = 42,
    similarity_k: int = 21,
    sketch_size: int = 2_000,
    max_train_similarity: float = 0.95,
    min_similarity_coverage: float = 0.85,
    similarity_bands: SimilarityBands = _DEFAULT_SIMILARITY_BANDS,
    similarity_table: Path | None = None,
) -> SplitPlan:
    """Create Test 2C and explicitly exclude candidates above the strict gate."""

    _validate_fraction(test_fraction)
    _validate_similarity_parameters(
        similarity_k=similarity_k,
        sketch_size=sketch_size,
        max_train_similarity=max_train_similarity,
        min_similarity_coverage=min_similarity_coverage,
        similarity_bands=similarity_bands,
    )

    groups = _prepare_groups(genomes)
    lookup = _group_lookup(groups)
    proposed = _stratified_partitions(
        groups,
        test_fraction=test_fraction,
        seed=seed,
        namespace="similarity:genome-proposal",
    )
    ordered = tuple(sorted(genomes, key=lambda genome: genome.genome_id))
    training = tuple(
        genome for genome in ordered if proposed[genome.genome_id] is SplitPartition.TRAIN
    )
    candidates = tuple(
        genome for genome in ordered if proposed[genome.genome_id] is SplitPartition.TEST
    )
    hits: dict[str, SimilarityHit]
    if similarity_table is None:
        hits = best_train_matches(
            training,
            candidates,
            k=similarity_k,
            sketch_size=sketch_size,
        )
        similarity_source = "built-in-minhash-mash"
        evidence_mode = "built-in-recomputable"
    else:
        hits = read_similarity_table(
            similarity_table,
            query_ids=(genome.genome_id for genome in candidates),
            reference_ids=(genome.genome_id for genome in training),
            max_train_similarity=max_train_similarity,
            min_similarity_coverage=min_similarity_coverage,
        )
        similarity_source = str(similarity_table)
        evidence_mode = "external-all-pairs"
    expected_queries = {genome.genome_id for genome in candidates}
    if set(hits) != expected_queries:
        missing = sorted(expected_queries - set(hits))
        extra = sorted(set(hits) - expected_queries)
        raise IntegrityError(
            "Similarity engine did not return exactly one best hit per candidate test genome; "
            f"missing={missing}, extra={extra}"
        )

    assignments: list[SplitAssignment] = []
    for genome in ordered:
        initial = proposed[genome.genome_id]
        group = lookup[genome.genome_id]
        if initial is SplitPartition.TRAIN:
            assignments.append(
                SplitAssignment(
                    genome_id=genome.genome_id,
                    label=genome.label,
                    partition=SplitPartition.TRAIN,
                    candidate_partition=SplitPartition.TRAIN,
                    reason="label_stratified_training_proposal",
                    group_id=group.group_id,
                    release_date=_release_date(genome),
                )
            )
            continue
        assignments.append(
            _similarity_candidate_assignment(
                genome,
                group_id=group.group_id,
                hit=hits[genome.genome_id],
                max_train_similarity=max_train_similarity,
                min_similarity_coverage=min_similarity_coverage,
                similarity_bands=similarity_bands,
            )
        )

    return SplitPlan(
        kind=SplitKind.SIMILARITY,
        assignments=tuple(assignments),
        seed=seed,
        parameters=_parameters(
            {
                "grouping": "canonical_topology_aware_genome_sha256_v2",
                "max_train_similarity": max_train_similarity,
                "min_similarity_coverage": min_similarity_coverage,
                "similarity_bands": similarity_bands,
                "similarity_k": similarity_k,
                "similarity_source": similarity_source,
                "similarity_evidence_mode": evidence_mode,
                "coverage_definition": (
                    "not_applicable_minhash"
                    if similarity_table is None
                    else "aligned_fraction_shorter"
                ),
                "similarity_table": similarity_table,
                "sketch_size": sketch_size,
                "strict_coverage_operator": ">=",
                "strict_identity_operator": ">",
                "test_fraction": test_fraction,
            }
        ),
    )


def _validate_partition(plan: SplitPlan) -> None:
    ids = [assignment.genome_id for assignment in plan.assignments]
    duplicates = sorted(genome_id for genome_id in set(ids) if ids.count(genome_id) > 1)
    if duplicates:
        raise IntegrityError("A genome has multiple split assignments: " + ", ".join(duplicates))
    if plan.train_ids & plan.test_ids or plan.train_ids & plan.excluded_ids:
        raise IntegrityError("Train, test, and excluded genome ID partitions are not disjoint")
    if plan.test_ids & plan.excluded_ids:
        raise IntegrityError("Train, test, and excluded genome ID partitions are not disjoint")

    by_group: dict[str, set[SplitPartition]] = defaultdict(set)
    for assignment in plan.assignments:
        if assignment.partition is not SplitPartition.EXCLUDED:
            by_group[assignment.group_id or assignment.genome_id].add(assignment.partition)
    leaking_groups = sorted(
        group_id
        for group_id, partitions in by_group.items()
        if SplitPartition.TRAIN in partitions and SplitPartition.TEST in partitions
    )
    if leaking_groups:
        raise IntegrityError(
            "Canonical genome content group occurs in both train and test: "
            + ", ".join(leaking_groups[:5])
        )

    for label in Label:
        train_count = sum(
            assignment.label is label and assignment.partition is SplitPartition.TRAIN
            for assignment in plan.assignments
        )
        test_count = sum(
            assignment.label is label and assignment.partition is SplitPartition.TEST
            for assignment in plan.assignments
        )
        if train_count == 0 or test_count == 0:
            raise IntegrityError(
                f"{plan.kind.value} split is not class-viable: {label.value} has "
                f"train={train_count}, test={test_count}"
            )


def _validate_temporal(plan: SplitPlan) -> None:
    cutoff = plan.parameters.get("release_date")
    if type(cutoff) is not date:
        raise IntegrityError("Temporal plan has no resolved release_date cutoff")
    effective: dict[str, date] = {}
    for assignment in plan.assignments:
        if assignment.release_date is not None:
            group_id = assignment.group_id or assignment.genome_id
            effective[group_id] = min(
                assignment.release_date,
                effective.get(group_id, assignment.release_date),
            )
    for assignment in plan.assignments:
        if assignment.release_date is None:
            if assignment.partition is not SplitPartition.EXCLUDED:
                raise IntegrityError(
                    f"Genome {assignment.genome_id!r} lacks release_date but is not excluded"
                )
            continue
        group_date = effective[assignment.group_id or assignment.genome_id]
        expected = SplitPartition.TRAIN if group_date <= cutoff else SplitPartition.TEST
        if assignment.partition is not expected:
            raise IntegrityError(
                f"Genome {assignment.genome_id!r} violates inclusive temporal cutoff "
                f"{cutoff.isoformat()}: group release_date={group_date.isoformat()}, "
                f"partition={assignment.partition.value}"
            )


def _validate_taxonomy(plan: SplitPlan) -> None:
    rank = plan.parameters.get("taxonomy_rank")
    holdouts = plan.parameters.get("holdout_taxa")
    if not isinstance(rank, str) or not isinstance(holdouts, tuple) or not holdouts:
        raise IntegrityError("Taxonomic plan lacks resolved rank and holdout_taxa")
    selected = {str(taxon).casefold() for taxon in holdouts}
    for assignment in plan.assignments:
        if assignment.label is Label.HOST:
            continue
        if assignment.taxon is None:
            if assignment.partition is not SplitPartition.EXCLUDED:
                raise IntegrityError(
                    f"Viral genome {assignment.genome_id!r} lacks {rank} but is not excluded"
                )
            continue
        is_holdout = assignment.taxon.casefold() in selected
        if is_holdout and assignment.partition is not SplitPartition.TEST:
            raise IntegrityError(
                f"Held-out {rank} {assignment.taxon!r} occurs outside the viral test partition"
            )
        if not is_holdout and assignment.partition is not SplitPartition.TRAIN:
            raise IntegrityError(
                f"Non-held-out viral {rank} {assignment.taxon!r} occurs outside training"
            )


def _validate_similarity(plan: SplitPlan) -> None:
    threshold = plan.parameters.get("max_train_similarity")
    min_coverage = plan.parameters.get("min_similarity_coverage")
    bands = plan.parameters.get("similarity_bands")
    if (
        not isinstance(threshold, (int, float))
        or isinstance(threshold, bool)
        or not isinstance(min_coverage, (int, float))
        or isinstance(min_coverage, bool)
        or not isinstance(bands, SimilarityBands)
    ):
        raise IntegrityError(
            "Similarity plan lacks resolved identity threshold, coverage threshold, "
            "or novelty bands"
        )
    train_ids = plan.train_ids
    for assignment in plan.assignments:
        if assignment.candidate_partition is None:
            raise IntegrityError(
                f"Similarity assignment {assignment.genome_id!r} lacks candidate_partition"
            )
        if assignment.candidate_partition is SplitPartition.TRAIN:
            if assignment.partition is not SplitPartition.TRAIN:
                raise IntegrityError("A proposed training genome changed partition")
            if any(
                value is not None
                for value in (
                    assignment.max_train_similarity,
                    assignment.nearest_train_genome_id,
                    assignment.similarity_method,
                    assignment.strict_gate_similarity,
                )
            ):
                raise IntegrityError("A proposed training genome carries test similarity evidence")
            continue
        if not assignment.similarity_method:
            raise IntegrityError(
                f"Candidate test genome {assignment.genome_id!r} lacks a similarity result"
            )
        expected_bin = bands.classify(assignment.max_train_similarity)
        if assignment.similarity_bin != expected_bin:
            raise IntegrityError(
                f"Genome {assignment.genome_id!r} has incorrect similarity bin "
                f"{assignment.similarity_bin!r}; expected {expected_bin!r}"
            )
        if (
            assignment.nearest_train_genome_id is not None
            and assignment.nearest_train_genome_id not in train_ids
        ):
            raise IntegrityError(
                f"Genome {assignment.genome_id!r} names a non-training nearest match"
            )
        if (assignment.max_train_similarity is None) != (
            assignment.nearest_train_genome_id is None
        ):
            raise IntegrityError(
                f"Genome {assignment.genome_id!r} has inconsistent maximum-hit provenance"
            )
        excluded_by_gate = assignment.strict_gate_similarity is not None
        maximum_itself_qualifies = (
            assignment.max_train_similarity is not None
            and assignment.max_train_similarity > threshold
            and (
                assignment.similarity_coverage is None
                or assignment.similarity_coverage >= min_coverage
            )
        )
        if maximum_itself_qualifies and not excluded_by_gate:
            raise IntegrityError(
                f"Genome {assignment.genome_id!r} violates the strict similarity gate: "
                "its qualifying maximum hit is missing strict-gate provenance"
            )
        if excluded_by_gate:
            if (
                assignment.strict_gate_similarity is None
                or assignment.strict_gate_similarity <= threshold
                or (
                    assignment.strict_gate_coverage is not None
                    and assignment.strict_gate_coverage < min_coverage
                )
                or assignment.strict_gate_train_genome_id not in train_ids
                or not assignment.strict_gate_method
            ):
                raise IntegrityError(
                    f"Genome {assignment.genome_id!r} has invalid strict-gate evidence"
                )
            if (
                assignment.max_train_similarity is not None
                and assignment.strict_gate_similarity > assignment.max_train_similarity
            ):
                raise IntegrityError(
                    f"Genome {assignment.genome_id!r} gate hit exceeds its recorded maximum"
                )
        expected_partition = SplitPartition.EXCLUDED if excluded_by_gate else SplitPartition.TEST
        if assignment.partition is not expected_partition:
            raise IntegrityError(
                f"Genome {assignment.genome_id!r} violates strict similarity gate: "
                f"identity > {threshold} with coverage >= {min_coverage} or absent"
            )


def build_split_plan(
    kind: SplitKind | str,
    genomes: Sequence[Genome],
    *,
    test_fraction: float = 0.20,
    seed: int = 42,
    missing_metadata: MissingMetadataPolicy | str = MissingMetadataPolicy.ERROR,
    temporal_cutoff: date | None = None,
    taxonomy_rank: str = "family",
    holdout_taxa: Iterable[str] = (),
    auto_holdout_count: int = 1,
    similarity_k: int = 21,
    sketch_size: int = 2_000,
    max_train_similarity: float = 0.95,
    min_similarity_coverage: float = 0.85,
    similarity_bands: SimilarityBands = _DEFAULT_SIMILARITY_BANDS,
    similarity_table: Path | None = None,
    historical_snapshot: str | None = None,
) -> SplitPlan:
    """Dispatch one source-genome protocol from resolved CLI/config values."""

    try:
        split_kind = SplitKind(kind)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"Unknown split kind {kind!r}") from exc
    if split_kind is SplitKind.RANDOM:
        raise ConfigurationError(
            "Random-fragment Test 2A is assigned after fragment generation and has no "
            "source-genome SplitPlan"
        )
    if split_kind is SplitKind.GENOME:
        return genome_holdout(genomes, test_fraction=test_fraction, seed=seed)
    if split_kind is SplitKind.TEMPORAL:
        return temporal_holdout(
            genomes,
            test_fraction=test_fraction,
            seed=seed,
            temporal_cutoff=temporal_cutoff,
            missing_metadata=missing_metadata,
            historical_snapshot=historical_snapshot,
        )
    if split_kind is SplitKind.TAXONOMY:
        return taxonomic_holdout(
            genomes,
            test_fraction=test_fraction,
            seed=seed,
            taxonomy_rank=taxonomy_rank,
            holdout_taxa=holdout_taxa,
            auto_holdout_count=auto_holdout_count,
            missing_metadata=missing_metadata,
        )
    return similarity_filtered_holdout(
        genomes,
        test_fraction=test_fraction,
        seed=seed,
        similarity_k=similarity_k,
        sketch_size=sketch_size,
        max_train_similarity=max_train_similarity,
        min_similarity_coverage=min_similarity_coverage,
        similarity_bands=similarity_bands,
        similarity_table=similarity_table,
    )


# Readable aliases for library users; the canonical names above mirror the CLI.
genome_level_split = genome_holdout
temporal_split = temporal_holdout
taxonomic_split = taxonomic_holdout
similarity_filtered_split = similarity_filtered_holdout


__all__ = [
    "FrozenParameters",
    "SplitAssignment",
    "SplitPartition",
    "SplitPlan",
    "build_split_plan",
    "genome_holdout",
    "genome_level_split",
    "similarity_filtered_holdout",
    "similarity_filtered_split",
    "taxonomic_holdout",
    "taxonomic_split",
    "temporal_holdout",
    "temporal_split",
]
