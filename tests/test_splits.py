from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
from datetime import date

import pytest

import chimera.splits as split_module
from chimera.config import MissingMetadataPolicy, SimilarityBands, SplitKind
from chimera.errors import ConfigurationError, InputError, IntegrityError
from chimera.models import Contig, Genome, GenomeMetadata, Label
from chimera.similarity import SimilarityHit
from chimera.splits import (
    FrozenParameters,
    SplitAssignment,
    SplitPartition,
    SplitPlan,
    build_split_plan,
    genome_holdout,
    similarity_filtered_holdout,
    taxonomic_holdout,
    temporal_holdout,
)


def _sequence(index: int) -> str:
    alphabet = "ACGT"
    encoded = "".join(alphabet[(index // (4**power)) % 4] for power in range(10))
    return ("ACGTTGCA" * 6) + encoded + "GATTACACCGTA"


def _metadata(
    *,
    release_date: date | None = None,
    family: str | None = None,
) -> GenomeMetadata:
    taxonomy = (("family", family),) if family is not None else ()
    return GenomeMetadata(release_date=release_date, taxonomy=taxonomy)


def _genome(
    genome_id: str,
    label: Label,
    index: int,
    *,
    sequence: str | None = None,
    release_date: date | None = None,
    family: str | None = None,
) -> Genome:
    return Genome(
        genome_id=genome_id,
        label=label,
        contigs=(Contig(f"contig-{genome_id}", sequence or _sequence(index)),),
        metadata=_metadata(release_date=release_date, family=family),
    )


def _balanced_genomes(per_class: int = 4) -> tuple[Genome, ...]:
    viruses = tuple(
        _genome(f"virus-{index}", Label.VIRUS, index) for index in range(1, per_class + 1)
    )
    hosts = tuple(
        _genome(f"host-{index}", Label.HOST, 100 + index) for index in range(1, per_class + 1)
    )
    return viruses + hosts


def _partition_by_id(plan: SplitPlan) -> dict[str, SplitPartition]:
    return {assignment.genome_id: assignment.partition for assignment in plan.assignments}


def _dated_genomes(genomes: tuple[Genome, ...]) -> tuple[Genome, ...]:
    return tuple(
        replace(
            genome,
            metadata=replace(genome.metadata, release_date=date(2020, 1, index)),
        )
        for index, genome in enumerate(genomes, start=1)
    )


def _assert_assignment_dates_match_references(
    plan: SplitPlan,
    genomes: tuple[Genome, ...],
) -> None:
    expected = {genome.genome_id: genome.metadata.release_date for genome in genomes}
    assert {
        assignment.genome_id: assignment.release_date for assignment in plan.assignments
    } == expected


def test_frozen_parameters_normalize_common_container_values() -> None:
    parameters = FrozenParameters.from_mapping(
        {
            "list": ["b", "a"],
            "set": {"b", "a"},
            "scalar": 3,
        }
    )

    assert tuple(parameters) == ("list", "scalar", "set")
    assert len(parameters) == 3
    assert parameters["list"] == ("b", "a")
    assert parameters["set"] == ("a", "b")
    with pytest.raises(KeyError):
        _ = parameters["missing"]


@pytest.mark.parametrize(
    ("entries", "error", "message"),
    [
        ((("", 1),), ValueError, "empty or duplicated"),
        ((("same", 1), ("same", 2)), ValueError, "empty or duplicated"),
        ((("unsupported", object()),), TypeError, "unsupported value type"),
    ],
)
def test_frozen_parameters_reject_invalid_entries(
    entries: tuple[tuple[str, object], ...],
    error: type[Exception],
    message: str,
) -> None:
    with pytest.raises(error, match=message):
        FrozenParameters(entries)  # type: ignore[arg-type]


def _base_assignment() -> dict[str, object]:
    return {
        "genome_id": "genome-a",
        "label": Label.VIRUS,
        "partition": SplitPartition.TRAIN,
        "reason": "curated assignment",
    }


@pytest.mark.parametrize(
    ("updates", "error", "message"),
    [
        ({"label": "bacterium"}, ValueError, "virus.*host"),
        ({"partition": "both"}, ValueError, "partition must"),
        ({"candidate_partition": "excluded"}, ValueError, "cannot be 'excluded'"),
        ({"candidate_partition": "both"}, ValueError, "candidate_partition must"),
        ({"genome_id": ""}, ValueError, "genome_id"),
        ({"reason": " "}, ValueError, "reason"),
        ({"reason": "two\nlines"}, ValueError, "one line"),
        ({"group_id": ""}, ValueError, "group_id"),
        ({"taxon": " "}, ValueError, "taxon"),
        ({"release_date": "2020-01-01"}, TypeError, "release_date"),
        ({"similarity_bin": "near"}, ValueError, "unknown similarity bin"),
        ({"max_train_similarity": float("nan")}, ValueError, "max_train_similarity"),
        ({"similarity_coverage": 1.1}, ValueError, "similarity_coverage"),
        ({"nearest_train_genome_id": ""}, ValueError, "nearest_train_genome_id"),
        ({"similarity_method": " "}, ValueError, "similarity_method"),
        (
            {"strict_gate_train_genome_id": "train-a"},
            ValueError,
            "requires strict_gate_similarity",
        ),
        (
            {"strict_gate_similarity": 0.99},
            ValueError,
            "requires a training genome ID and method",
        ),
    ],
)
def test_split_assignment_rejects_invalid_serializable_state(
    updates: dict[str, object],
    error: type[Exception],
    message: str,
) -> None:
    with pytest.raises(error, match=message):
        SplitAssignment(**(_base_assignment() | updates))  # type: ignore[arg-type]


def test_split_assignment_normalizes_defaults_and_text() -> None:
    assignment = SplitAssignment(**(_base_assignment() | {"reason": "  curated  ", "taxon": " X "}))  # type: ignore[arg-type]

    assert assignment.group_id == "genome-a"
    assert assignment.reason == "curated"
    assert assignment.taxon == "X"


def test_genome_holdout_is_semantically_seeded_and_class_stratified() -> None:
    genomes = _balanced_genomes()

    forward = genome_holdout(genomes, test_fraction=0.5, seed=8128)
    reversed_input = genome_holdout(tuple(reversed(genomes)), test_fraction=0.5, seed=8128)

    assert forward == reversed_input
    assert forward.kind is SplitKind.GENOME
    assert {assignment.genome_id for assignment in forward.assignments} == {
        genome.genome_id for genome in genomes
    }
    for label in Label:
        assert any(item.label is label for item in forward.train)
        assert any(item.label is label for item in forward.test)
    assert forward.train_ids.isdisjoint(forward.test_ids)
    assert forward.parameters["test_fraction"] == 0.5


def test_non_temporal_assignments_preserve_effective_reference_dates() -> None:
    genome_inputs = _dated_genomes(_balanced_genomes())
    taxonomy_inputs = _dated_genomes(_taxonomy_genomes(include_missing=False))

    genome_plan = genome_holdout(genome_inputs, test_fraction=0.5, seed=17)
    taxonomy_plan = taxonomic_holdout(
        taxonomy_inputs,
        holdout_taxa=("Alpha",),
        test_fraction=0.5,
        seed=19,
    )
    similarity_plan = similarity_filtered_holdout(
        genome_inputs,
        test_fraction=0.5,
        seed=23,
    )

    _assert_assignment_dates_match_references(genome_plan, genome_inputs)
    _assert_assignment_dates_match_references(taxonomy_plan, taxonomy_inputs)
    _assert_assignment_dates_match_references(similarity_plan, genome_inputs)
    assert {assignment.candidate_partition for assignment in similarity_plan.assignments} == {
        SplitPartition.TRAIN,
        SplitPartition.TEST,
    }


def test_split_records_and_parameters_are_immutable() -> None:
    plan = genome_holdout(_balanced_genomes(), seed=7)

    with pytest.raises(FrozenInstanceError):
        plan.assignments[0].reason = "changed"  # type: ignore[misc]
    with pytest.raises(FrozenInstanceError):
        plan.seed = 99  # type: ignore[misc]
    with pytest.raises(TypeError):
        plan.parameters["test_fraction"] = 0.9  # type: ignore[index]


def test_split_plan_constructor_rejects_invalid_core_fields() -> None:
    assignment = SplitAssignment(**_base_assignment())  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="supported SplitKind"):
        SplitPlan("unknown", (assignment,), 1)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="seed must be an integer"):
        SplitPlan(SplitKind.GENOME, (assignment,), True)
    with pytest.raises(TypeError, match="iterable of SplitAssignment"):
        SplitPlan(SplitKind.GENOME, "not-assignments", 1)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="contain SplitAssignment"):
        SplitPlan(SplitKind.GENOME, (object(),), 1)  # type: ignore[arg-type]


def test_split_plan_parameter_coercion_and_unknown_assignment_lookup() -> None:
    plan = genome_holdout(_balanced_genomes(), test_fraction=0.5, seed=7)
    replaced = replace(plan, parameters={"custom": ["value"]})  # type: ignore[arg-type]

    assert isinstance(replaced.parameters, FrozenParameters)
    assert replaced.parameters["custom"] == ("value",)
    with pytest.raises(KeyError, match="unknown-genome"):
        replaced.assignment_for("unknown-genome")
    with pytest.raises(TypeError, match="parameters must be a mapping"):
        replace(plan, parameters=object())  # type: ignore[arg-type]


@pytest.mark.parametrize("test_fraction", [0.0, 1.0, float("nan"), float("inf")])
def test_genome_holdout_rejects_invalid_test_fraction(test_fraction: float) -> None:
    with pytest.raises(ConfigurationError, match="finite and strictly between"):
        genome_holdout(_balanced_genomes(), test_fraction=test_fraction)


@pytest.mark.parametrize(
    ("genomes", "error", "message"),
    [
        ("not-genomes", TypeError, "sequence of Genome"),
        ((), ConfigurationError, "At least one source genome"),
        ((object(),), TypeError, "only Genome"),
    ],
)
def test_genome_holdout_rejects_malformed_source_collections(
    genomes: object,
    error: type[Exception],
    message: str,
) -> None:
    with pytest.raises(error, match=message):
        genome_holdout(genomes)  # type: ignore[arg-type]


def test_genome_holdout_rejects_duplicate_genome_identifiers() -> None:
    genomes = list(_balanced_genomes())
    genomes.append(_genome(genomes[0].genome_id, Label.VIRUS, 999))

    with pytest.raises(InputError, match="Duplicate genome_id"):
        genome_holdout(genomes)


def test_exact_duplicate_genomes_cannot_cross_the_partition() -> None:
    genomes = list(_balanced_genomes())
    duplicated_sequence = genomes[0].contigs[0].sequence
    genomes.append(_genome("virus-duplicate", Label.VIRUS, 999, sequence=duplicated_sequence))

    plan = genome_holdout(genomes, test_fraction=0.5, seed=19)
    original = plan.assignment_for(genomes[0].genome_id)
    duplicate = plan.assignment_for("virus-duplicate")

    assert original.group_id == duplicate.group_id
    assert original.partition is duplicate.partition


def test_validator_rejects_a_duplicate_content_leakage_trap() -> None:
    genomes = list(_balanced_genomes())
    genomes.append(
        _genome(
            "virus-duplicate",
            Label.VIRUS,
            999,
            sequence=genomes[0].contigs[0].sequence,
        )
    )
    plan = genome_holdout(genomes, test_fraction=0.5, seed=29)
    first = plan.assignment_for("virus-1")
    duplicate = plan.assignment_for("virus-duplicate")
    opposite = (
        SplitPartition.TEST if first.partition is SplitPartition.TRAIN else SplitPartition.TRAIN
    )
    corrupted = replace(duplicate, partition=opposite)
    assignments = tuple(
        corrupted if item.genome_id == corrupted.genome_id else item for item in plan.assignments
    )

    with pytest.raises(IntegrityError, match="content group"):
        replace(plan, assignments=assignments)


def test_conflicting_labels_for_identical_content_are_quarantined() -> None:
    sequence = _sequence(1)
    genomes = list(_balanced_genomes())
    genomes[4] = _genome("host-conflict", Label.HOST, 101, sequence=sequence)

    with pytest.raises(InputError, match="conflicting virus/host labels"):
        genome_holdout(genomes)


def test_genome_holdout_requires_two_independent_groups_per_class() -> None:
    genomes = (
        _genome("virus-a", Label.VIRUS, 1),
        _genome("host-a", Label.HOST, 101),
        _genome("host-b", Label.HOST, 102),
    )

    with pytest.raises(ConfigurationError, match="at least two independent"):
        genome_holdout(genomes)


def _temporal_genomes(*, include_missing: bool = True) -> tuple[Genome, ...]:
    days = (date(2020, 1, 1), date(2021, 1, 1), date(2022, 1, 1), date(2023, 1, 1))
    result = [
        _genome(f"virus-{index}", Label.VIRUS, index, release_date=value)
        for index, value in enumerate(days, start=1)
    ]
    result.extend(
        _genome(f"host-{index}", Label.HOST, 100 + index, release_date=value)
        for index, value in enumerate(days, start=1)
    )
    # A later accession with already-public sequence content belongs in train.
    result.append(
        _genome(
            "virus-late-duplicate",
            Label.VIRUS,
            999,
            sequence=result[0].contigs[0].sequence,
            release_date=date(2024, 1, 1),
        )
    )
    if include_missing:
        result.append(_genome("virus-undated", Label.VIRUS, 77))
    return tuple(result)


def test_temporal_cutoff_is_inclusive_and_missing_dates_never_train() -> None:
    cutoff = date(2021, 1, 1)
    genomes = _temporal_genomes()

    plan = temporal_holdout(
        genomes,
        temporal_cutoff=cutoff,
        missing_metadata=MissingMetadataPolicy.EXCLUDE,
    )

    assert plan.parameters["release_date"] == cutoff
    assert plan.parameters["temporal_semantics"] == "release-date-filtered retrospective"
    assert plan.assignment_for("virus-2").partition is SplitPartition.TRAIN
    assert plan.assignment_for("host-2").partition is SplitPartition.TRAIN
    assert plan.assignment_for("virus-3").partition is SplitPartition.TEST
    assert plan.assignment_for("virus-undated").partition is SplitPartition.EXCLUDED
    assert plan.assignment_for("virus-undated").reason == "missing_release_date"
    assert plan.assignment_for("virus-late-duplicate").partition is SplitPartition.TRAIN
    assert {item.genome_id for item in plan.assignments} == {genome.genome_id for genome in genomes}


def test_temporal_auto_cutoff_is_viable_and_input_order_independent() -> None:
    genomes = _temporal_genomes(include_missing=False)

    first = temporal_holdout(genomes, test_fraction=0.5, seed=11)
    second = temporal_holdout(tuple(reversed(genomes)), test_fraction=0.5, seed=11)

    assert first == second
    assert type(first.parameters["release_date"]) is date
    assert first.parameters["cutoff_selection"] == "automatic"
    for label in Label:
        assert any(item.label is label for item in first.train)
        assert any(item.label is label for item in first.test)


def test_temporal_missing_metadata_error_is_actionable() -> None:
    with pytest.raises(ConfigurationError, match="lack release_date"):
        temporal_holdout(_temporal_genomes(), missing_metadata=MissingMetadataPolicy.ERROR)


def test_temporal_parameter_validation_and_snapshot_semantics() -> None:
    genomes = _temporal_genomes(include_missing=False)
    cutoff = date(2021, 1, 1)

    with pytest.raises(ConfigurationError, match="disagree"):
        temporal_holdout(
            genomes,
            release_date=cutoff,
            temporal_cutoff=date(2022, 1, 1),
        )
    with pytest.raises(ConfigurationError, match=r"datetime\.date"):
        temporal_holdout(genomes, release_date="2021-01-01")  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError, match="historical_snapshot cannot be blank"):
        temporal_holdout(genomes, temporal_cutoff=cutoff, historical_snapshot="  ")
    with pytest.raises(ConfigurationError, match="missing_metadata"):
        temporal_holdout(genomes, missing_metadata="drop")

    plan = temporal_holdout(
        genomes,
        release_date=cutoff,
        temporal_cutoff=cutoff,
        historical_snapshot="NCBI Virus 2021-01-01 snapshot",
    )
    assert plan.parameters["temporal_semantics"] == "historical-snapshot prospective"


def test_temporal_auto_cutoff_rejects_nonviable_or_fully_undated_catalogs() -> None:
    same_date = date(2020, 1, 1)
    nonviable = tuple(
        _genome(f"virus-{index}", Label.VIRUS, index, release_date=same_date) for index in range(2)
    ) + tuple(
        _genome(f"host-{index}", Label.HOST, 100 + index, release_date=same_date)
        for index in range(2)
    )
    with pytest.raises(ConfigurationError, match="No release-date cutoff"):
        temporal_holdout(nonviable)

    with pytest.raises(ConfigurationError, match="No genomes have release_date"):
        temporal_holdout(
            _balanced_genomes(),
            missing_metadata=MissingMetadataPolicy.EXCLUDE,
        )


def test_temporal_validator_catches_a_cutoff_leakage_trap() -> None:
    plan = temporal_holdout(
        _temporal_genomes(include_missing=False),
        temporal_cutoff=date(2021, 1, 1),
    )
    leaked = replace(plan.assignment_for("virus-3"), partition=SplitPartition.TRAIN)
    assignments = tuple(
        leaked if item.genome_id == leaked.genome_id else item for item in plan.assignments
    )

    with pytest.raises(IntegrityError, match="temporal cutoff"):
        replace(plan, assignments=assignments)


def test_temporal_validator_requires_cutoff_and_exclusion_of_undated_sources() -> None:
    plan = temporal_holdout(
        _temporal_genomes(),
        temporal_cutoff=date(2021, 1, 1),
        missing_metadata=MissingMetadataPolicy.EXCLUDE,
    )
    with pytest.raises(IntegrityError, match="no resolved release_date cutoff"):
        replace(plan, parameters={})  # type: ignore[arg-type]

    undated = replace(
        plan.assignment_for("virus-undated"),
        partition=SplitPartition.TRAIN,
    )
    assignments = tuple(
        undated if item.genome_id == undated.genome_id else item for item in plan.assignments
    )
    with pytest.raises(IntegrityError, match="lacks release_date but is not excluded"):
        replace(plan, assignments=assignments)


def _taxonomy_genomes(*, include_missing: bool = True) -> tuple[Genome, ...]:
    viruses = (
        _genome("virus-alpha-1", Label.VIRUS, 1, family="Alpha"),
        _genome("virus-alpha-2", Label.VIRUS, 2, family="Alpha"),
        _genome("virus-beta", Label.VIRUS, 3, family="Beta"),
        _genome("virus-gamma", Label.VIRUS, 4, family="Gamma"),
    )
    hosts = tuple(_genome(f"host-{index}", Label.HOST, 100 + index) for index in range(1, 5))
    missing = (_genome("virus-unclassified", Label.VIRUS, 5),) if include_missing else ()
    return viruses + hosts + missing


def test_explicit_taxonomic_holdout_excludes_every_selected_viral_taxon() -> None:
    genomes = _taxonomy_genomes()
    plan = taxonomic_holdout(
        genomes,
        taxonomy_rank="family",
        holdout_taxa=("alpha",),
        test_fraction=0.5,
        missing_metadata=MissingMetadataPolicy.EXCLUDE,
    )

    assert plan.parameters["holdout_taxa"] == ("Alpha",)
    assert plan.assignment_for("virus-unclassified").partition is SplitPartition.EXCLUDED
    assert all(
        assignment.taxon != "Alpha" for assignment in plan.train if assignment.label is Label.VIRUS
    )
    assert {
        assignment.genome_id for assignment in plan.test if assignment.label is Label.VIRUS
    } == {"virus-alpha-1", "virus-alpha-2"}
    assert any(item.label is Label.HOST for item in plan.train)
    assert any(item.label is Label.HOST for item in plan.test)
    assert {item.genome_id for item in plan.assignments} == {genome.genome_id for genome in genomes}


def test_auto_taxon_selection_is_semantically_seeded() -> None:
    genomes = _taxonomy_genomes(include_missing=False)

    first = taxonomic_holdout(genomes, seed=123, auto_holdout_count=1)
    second = taxonomic_holdout(tuple(reversed(genomes)), seed=123, auto_holdout_count=1)

    assert first == second
    assert first.parameters["holdout_selection"] == "automatic"
    assert len(first.parameters["holdout_taxa"]) == 1  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"taxonomy_rank": ""}, "taxonomy_rank"),
        ({"taxonomy_rank": "viral family"}, "taxonomy_rank"),
        ({"auto_holdout_count": 0}, "auto_holdout_count"),
        ({"auto_holdout_count": True}, "auto_holdout_count"),
        ({"holdout_taxa": ("",)}, "blank names"),
    ],
)
def test_taxonomic_parameter_validation_is_actionable(
    updates: dict[str, object],
    message: str,
) -> None:
    with pytest.raises(ConfigurationError, match=message):
        taxonomic_holdout(_taxonomy_genomes(include_missing=False), **updates)  # type: ignore[arg-type]


def test_taxonomic_holdout_normalizes_duplicate_requested_names() -> None:
    plan = taxonomic_holdout(
        _taxonomy_genomes(include_missing=False),
        holdout_taxa=(" alpha ", "ALPHA"),
        test_fraction=0.5,
    )
    assert plan.parameters["holdout_taxa"] == ("Alpha",)


def test_taxonomic_holdout_rejects_missing_or_insufficient_taxonomy() -> None:
    with pytest.raises(ConfigurationError, match="lack taxonomy rank"):
        taxonomic_holdout(_taxonomy_genomes())

    one_taxon = (
        _genome("virus-a", Label.VIRUS, 1, family="Alpha"),
        _genome("virus-b", Label.VIRUS, 2, family="Alpha"),
        _genome("host-a", Label.HOST, 101),
        _genome("host-b", Label.HOST, 102),
    )
    with pytest.raises(ConfigurationError, match="at least two represented viral taxa"):
        taxonomic_holdout(one_taxon)


def test_taxonomic_holdout_must_leave_a_training_taxon() -> None:
    genomes = _taxonomy_genomes(include_missing=False)
    with pytest.raises(ConfigurationError, match="leave at least one represented viral taxon"):
        taxonomic_holdout(genomes, holdout_taxa=("Alpha", "Beta", "Gamma"))
    with pytest.raises(ConfigurationError, match="leave at least one represented viral taxon"):
        taxonomic_holdout(genomes, auto_holdout_count=3)


def test_taxonomy_validator_catches_a_heldout_taxon_in_training() -> None:
    plan = taxonomic_holdout(
        _taxonomy_genomes(include_missing=False),
        holdout_taxa=("Alpha",),
        test_fraction=0.5,
    )
    leaked = replace(plan.assignment_for("virus-alpha-1"), partition=SplitPartition.TRAIN)
    assignments = tuple(
        leaked if item.genome_id == leaked.genome_id else item for item in plan.assignments
    )

    with pytest.raises(IntegrityError, match="Held-out family"):
        replace(plan, assignments=assignments)


def test_taxonomy_validator_rejects_missing_parameters_and_invalid_viral_states() -> None:
    plan = taxonomic_holdout(
        _taxonomy_genomes(),
        holdout_taxa=("Alpha",),
        test_fraction=0.5,
        missing_metadata=MissingMetadataPolicy.EXCLUDE,
    )
    with pytest.raises(IntegrityError, match="lacks resolved rank"):
        replace(plan, parameters={})  # type: ignore[arg-type]

    unclassified = replace(
        plan.assignment_for("virus-unclassified"),
        partition=SplitPartition.TRAIN,
    )
    assignments = tuple(
        unclassified if item.genome_id == unclassified.genome_id else item
        for item in plan.assignments
    )
    with pytest.raises(IntegrityError, match="lacks family but is not excluded"):
        replace(plan, assignments=assignments)

    non_holdout = replace(
        plan.assignment_for("virus-beta"),
        partition=SplitPartition.TEST,
    )
    assignments = tuple(
        non_holdout if item.genome_id == non_holdout.genome_id else item
        for item in plan.assignments
    )
    with pytest.raises(IntegrityError, match="Non-held-out viral family"):
        replace(plan, assignments=assignments)


def test_taxonomic_holdout_rejects_unknown_requested_taxon() -> None:
    with pytest.raises(ConfigurationError, match="absent"):
        taxonomic_holdout(
            _taxonomy_genomes(include_missing=False),
            holdout_taxa=("NotARealFamily",),
        )


def test_taxonomic_holdout_rejects_conflicting_taxonomy_for_duplicate_content() -> None:
    sequence = _sequence(1)
    genomes = list(_taxonomy_genomes(include_missing=False))
    genomes[1] = _genome(
        "virus-conflict",
        Label.VIRUS,
        2,
        sequence=sequence,
        family="Beta",
    )

    with pytest.raises(InputError, match="conflicting family taxonomy"):
        taxonomic_holdout(genomes, holdout_taxa=("Alpha",))


def test_similarity_split_audits_every_candidate_and_never_silently_drops(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    genomes = _balanced_genomes()
    observed: dict[str, tuple[str, ...]] = {}

    def fake_best_matches(
        train: tuple[Genome, ...],
        test: tuple[Genome, ...],
        *,
        k: int,
        sketch_size: int,
    ) -> dict[str, SimilarityHit]:
        observed["train"] = tuple(genome.genome_id for genome in train)
        observed["test"] = tuple(genome.genome_id for genome in test)
        assert k == 21
        assert sketch_size == 2_000
        reference = train[0].genome_id
        virus_candidates = sorted(
            genome.genome_id for genome in test if genome.label is Label.VIRUS
        )
        first_virus = virus_candidates[0]
        low_coverage_virus = virus_candidates[-1]
        first_host = min(genome.genome_id for genome in test if genome.label is Label.HOST)
        result: dict[str, SimilarityHit] = {}
        for genome in test:
            similarity = 0.40
            coverage = 0.90
            if genome.genome_id == first_virus:
                similarity = 0.99
            elif genome.genome_id == low_coverage_virus:
                similarity = 0.99
                coverage = 0.50
            elif genome.genome_id == first_host:
                similarity = 0.95  # Equality is retained: the strict operator is >.
            result[genome.genome_id] = SimilarityHit(
                query_genome_id=genome.genome_id,
                reference_genome_id=reference,
                similarity=similarity,
                coverage=coverage,
                method="test-alignment",
            )
        return result

    monkeypatch.setattr(split_module, "best_train_matches", fake_best_matches)
    plan = similarity_filtered_holdout(
        genomes,
        test_fraction=0.5,
        seed=73,
        max_train_similarity=0.95,
        similarity_bands=SimilarityBands(high=0.9, moderate=0.7, low=0.3),
    )

    assert set(observed["train"]) == plan.train_ids
    candidate_test_ids = {
        item.genome_id
        for item in plan.assignments
        if item.candidate_partition is SplitPartition.TEST
    }
    assert set(observed["test"]) == candidate_test_ids
    assert {item.genome_id for item in plan.assignments} == {genome.genome_id for genome in genomes}
    assert len(plan.excluded) == 1
    excluded = plan.excluded[0]
    assert excluded.candidate_partition is SplitPartition.TEST
    assert excluded.reason == "similarity_above_strict_identity_and_coverage_gate"
    assert excluded.similarity_bin == "high_similarity"
    assert excluded.max_train_similarity == 0.99
    assert all(
        item.max_train_similarity is None
        or item.max_train_similarity <= 0.95
        or (item.similarity_coverage is not None and item.similarity_coverage < 0.85)
        for item in plan.test
    )
    assert any(item.max_train_similarity == 0.95 for item in plan.test)
    low_coverage = next(
        item
        for item in plan.test
        if item.max_train_similarity == 0.99 and item.similarity_coverage == 0.50
    )
    assert low_coverage.reason == "similarity_above_identity_below_minimum_coverage"
    assert plan.parameters["min_similarity_coverage"] == 0.85
    assert all(item.nearest_train_genome_id in plan.train_ids for item in plan.test + plan.excluded)

    retained = next(item for item in plan.test if item.max_train_similarity == 0.40)
    leaked = replace(
        retained,
        max_train_similarity=0.96,
        similarity_bin="high_similarity",
    )
    assignments = tuple(
        leaked if item.genome_id == leaked.genome_id else item for item in plan.assignments
    )
    with pytest.raises(IntegrityError, match="strict similarity gate"):
        replace(plan, assignments=assignments)

    now_sufficient = replace(low_coverage, similarity_coverage=0.85)
    coverage_leak = tuple(
        now_sufficient if item.genome_id == now_sufficient.genome_id else item
        for item in plan.assignments
    )
    with pytest.raises(IntegrityError, match="strict similarity gate"):
        replace(plan, assignments=coverage_leak)


def test_similarity_split_rejects_an_incomplete_candidate_audit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def incomplete_results(
        train: tuple[Genome, ...],
        test: tuple[Genome, ...],
        *,
        k: int,
        sketch_size: int,
    ) -> dict[str, SimilarityHit]:
        del train, k, sketch_size
        genome = test[0]
        return {
            genome.genome_id: SimilarityHit(
                query_genome_id=genome.genome_id,
                reference_genome_id=None,
                similarity=None,
                coverage=None,
                method="incomplete-test-engine",
            )
        }

    monkeypatch.setattr(split_module, "best_train_matches", incomplete_results)

    with pytest.raises(IntegrityError, match="exactly one best hit"):
        similarity_filtered_holdout(_balanced_genomes(), test_fraction=0.5)


@pytest.mark.parametrize(
    ("updates", "error", "message"),
    [
        ({"similarity_k": True}, ConfigurationError, "similarity_k"),
        ({"similarity_k": 0}, ConfigurationError, "similarity_k"),
        ({"sketch_size": True}, ConfigurationError, "sketch_size"),
        ({"sketch_size": 0}, ConfigurationError, "sketch_size"),
        ({"max_train_similarity": float("nan")}, ConfigurationError, "max_train_similarity"),
        ({"max_train_similarity": 1.1}, ConfigurationError, "max_train_similarity"),
        ({"min_similarity_coverage": float("inf")}, ConfigurationError, "min_similarity_coverage"),
        ({"min_similarity_coverage": -0.1}, ConfigurationError, "min_similarity_coverage"),
        ({"similarity_bands": object()}, TypeError, "SimilarityBands"),
        (
            {
                "max_train_similarity": 0.8,
                "similarity_bands": SimilarityBands(high=0.9, moderate=0.7, low=0.3),
            },
            ConfigurationError,
            "high novelty-band",
        ),
    ],
)
def test_similarity_parameter_validation_is_actionable(
    updates: dict[str, object],
    error: type[Exception],
    message: str,
) -> None:
    with pytest.raises(error, match=message):
        similarity_filtered_holdout(_balanced_genomes(), **updates)  # type: ignore[arg-type]


def test_similarity_external_table_state_is_recorded(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    table = tmp_path / "external.tsv"

    def fake_external_table(
        path,
        *,
        query_ids,
        reference_ids,
        max_train_similarity,
        min_similarity_coverage,
    ):
        queries = tuple(query_ids)
        references = tuple(reference_ids)
        assert path == table
        assert max_train_similarity == 0.95
        assert min_similarity_coverage == 0.85
        return {
            query: SimilarityHit(
                query_genome_id=query,
                reference_genome_id=references[0],
                similarity=0.4,
                coverage=0.9,
                method="skani-0.3.1",
            )
            for query in queries
        }

    monkeypatch.setattr(split_module, "read_similarity_table", fake_external_table)
    plan = similarity_filtered_holdout(
        _balanced_genomes(),
        test_fraction=0.5,
        similarity_table=table,
    )

    assert plan.parameters["similarity_source"] == str(table)
    assert plan.parameters["similarity_evidence_mode"] == "external-all-pairs"
    assert plan.parameters["coverage_definition"] == "aligned_fraction_shorter"


def test_similarity_rejects_hit_with_mismatched_query(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def mismatched_hits(train, test, *, k, sketch_size):
        del k, sketch_size
        reference = train[0].genome_id
        return {
            genome.genome_id: SimilarityHit(
                query_genome_id="wrong-query",
                reference_genome_id=reference,
                similarity=0.4,
                coverage=None,
                method="test-minhash-1",
            )
            for genome in test
        }

    monkeypatch.setattr(split_module, "best_train_matches", mismatched_hits)
    with pytest.raises(IntegrityError, match="names query"):
        similarity_filtered_holdout(_balanced_genomes(), test_fraction=0.5)


def test_partition_validator_rejects_duplicate_ids_and_nonviable_classes() -> None:
    plan = genome_holdout(_balanced_genomes(), test_fraction=0.5, seed=31)
    first, second, *remaining = plan.assignments
    duplicate = replace(second, genome_id=first.genome_id)
    with pytest.raises(IntegrityError, match="multiple split assignments"):
        replace(plan, assignments=(first, duplicate, *remaining))

    no_virus_test = tuple(
        replace(assignment, partition=SplitPartition.EXCLUDED)
        if assignment.label is Label.VIRUS and assignment.partition is SplitPartition.TEST
        else assignment
        for assignment in plan.assignments
    )
    with pytest.raises(IntegrityError, match="not class-viable: virus"):
        replace(plan, assignments=no_virus_test)


def test_similarity_validator_rejects_corrupted_protocol_states() -> None:
    plan = similarity_filtered_holdout(_balanced_genomes(), test_fraction=0.5, seed=37)
    training = plan.train[0]
    candidate = next(
        assignment
        for assignment in plan.assignments
        if assignment.candidate_partition is SplitPartition.TEST
    )
    training_id = training.genome_id

    def assert_rejected(changed: SplitAssignment, message: str) -> None:
        assignments = tuple(
            changed if item.genome_id == changed.genome_id else item for item in plan.assignments
        )
        with pytest.raises(IntegrityError, match=message):
            replace(plan, assignments=assignments)

    with pytest.raises(IntegrityError, match="lacks resolved identity threshold"):
        replace(plan, parameters={})  # type: ignore[arg-type]

    assert_rejected(replace(candidate, candidate_partition=None), "lacks candidate_partition")
    replacement_train_id = next(
        assignment.genome_id for assignment in plan.train if assignment.genome_id != training_id
    )
    changed_training_state = tuple(
        replace(item, partition=SplitPartition.TEST)
        if item.genome_id == training_id
        else replace(
            item,
            nearest_train_genome_id=(
                replacement_train_id
                if item.nearest_train_genome_id == training_id
                else item.nearest_train_genome_id
            ),
            strict_gate_train_genome_id=(
                replacement_train_id
                if item.strict_gate_train_genome_id == training_id
                else item.strict_gate_train_genome_id
            ),
        )
        for item in plan.assignments
    )
    with pytest.raises(IntegrityError, match="proposed training genome changed partition"):
        replace(plan, assignments=changed_training_state)
    assert_rejected(
        replace(training, similarity_method="unexpected-1"),
        "training genome carries test similarity evidence",
    )
    assert_rejected(
        replace(candidate, similarity_method=None),
        "lacks a similarity result",
    )
    wrong_bin = (
        "low_similarity" if candidate.similarity_bin != "low_similarity" else "high_similarity"
    )
    assert_rejected(replace(candidate, similarity_bin=wrong_bin), "incorrect similarity bin")
    assert_rejected(
        replace(candidate, nearest_train_genome_id="not-a-training-genome"),
        "non-training nearest match",
    )
    assert_rejected(
        replace(
            candidate,
            max_train_similarity=0.4,
            nearest_train_genome_id=None,
            similarity_bin="low_similarity",
        ),
        "inconsistent maximum-hit provenance",
    )

    invalid_gate = replace(
        candidate,
        partition=SplitPartition.EXCLUDED,
        max_train_similarity=0.95,
        nearest_train_genome_id=training_id,
        similarity_coverage=0.9,
        similarity_method="alignment-1",
        similarity_bin="high_similarity",
        strict_gate_train_genome_id=training_id,
        strict_gate_similarity=0.95,
        strict_gate_coverage=0.9,
        strict_gate_method="alignment-1",
    )
    assert_rejected(invalid_gate, "invalid strict-gate evidence")

    gate_above_maximum = replace(
        invalid_gate,
        max_train_similarity=0.98,
        strict_gate_similarity=0.99,
    )
    assert_rejected(gate_above_maximum, "gate hit exceeds its recorded maximum")

    retained_despite_gate = replace(
        invalid_gate,
        partition=SplitPartition.TEST,
        max_train_similarity=0.98,
        strict_gate_similarity=0.98,
    )
    assert_rejected(retained_despite_gate, "violates strict similarity gate")


def test_build_split_plan_dispatch_and_random_fragment_boundary() -> None:
    plan = build_split_plan(
        SplitKind.GENOME,
        _balanced_genomes(),
        test_fraction=0.5,
        seed=3,
    )
    assert plan.kind is SplitKind.GENOME

    with pytest.raises(ConfigurationError, match="Random-fragment"):
        build_split_plan(SplitKind.RANDOM, _balanced_genomes())


def test_build_split_plan_dispatches_all_source_protocols() -> None:
    temporal = build_split_plan(
        "temporal",
        _temporal_genomes(include_missing=False),
        temporal_cutoff=date(2021, 1, 1),
    )
    taxonomy = build_split_plan(
        "taxonomy",
        _taxonomy_genomes(include_missing=False),
        holdout_taxa=("Alpha",),
        test_fraction=0.5,
    )
    similarity = build_split_plan(
        "similarity",
        _balanced_genomes(),
        test_fraction=0.5,
    )

    assert temporal.kind is SplitKind.TEMPORAL
    assert taxonomy.kind is SplitKind.TAXONOMY
    assert similarity.kind is SplitKind.SIMILARITY
    with pytest.raises(ConfigurationError, match="Unknown split kind"):
        build_split_plan("unknown", _balanced_genomes())
