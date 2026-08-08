"""Tests for deterministic, leakage-safe fragment generation."""

from __future__ import annotations

from collections import Counter, defaultdict

import pytest

import chimera.fragments as fragment_module
from chimera.errors import InputError, IntegrityError
from chimera.fragments import derive_seed, generate_fragments, split_fragments_random
from chimera.models import Contig, Fragment, Genome, Label, reverse_complement


def _genome(genome_id: str, label: Label, *sequences: str) -> Genome:
    return Genome(
        genome_id=genome_id,
        label=label,
        contigs=tuple(
            Contig(sequence_id=f"{genome_id}.contig{index}", sequence=sequence)
            for index, sequence in enumerate(sequences, start=1)
        ),
    )


def _source_sequence(fragment: Fragment, genomes: tuple[Genome, ...]) -> str:
    contigs = {
        (genome.genome_id, contig.sequence_id): contig.sequence
        for genome in genomes
        for contig in genome.contigs
    }
    return contigs[(fragment.genome_id, fragment.sequence_id)][fragment.start : fragment.end]


def test_derive_seed_is_stable_typed_and_semantic() -> None:
    assert derive_seed(42, "coordinates", "genome-1", 500) == derive_seed(
        42, "coordinates", "genome-1", 500
    )
    assert derive_seed(42, "ab", "c") != derive_seed(42, "a", "bc")
    assert derive_seed(42, "1") != derive_seed(42, 1)
    assert derive_seed(42, "coordinates") != derive_seed(42, "strands")

    with pytest.raises(TypeError, match="seed must be an integer"):
        derive_seed(True, "coordinates")
    with pytest.raises(TypeError, match="strings or integers"):
        derive_seed(42, 1.5)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="not bool"):
        derive_seed(42, True)


@pytest.mark.parametrize(
    ("genomes", "message"),
    [
        ("not-genomes", "iterable of Genome"),
        (42, "iterable of Genome"),
        ((), "at least one genome"),
        ((object(),), "only Genome"),
    ],
)
def test_generation_rejects_malformed_genome_collections(
    genomes: object,
    message: str,
) -> None:
    with pytest.raises(InputError, match=message):
        generate_fragments(
            genomes,  # type: ignore[arg-type]
            fragment_lengths=(4,),
            fragments_per_genome=2,
            seed=1,
        )


def test_generation_rejects_duplicate_genome_ids() -> None:
    genomes = (
        _genome("duplicate", Label.VIRUS, "ACGT" * 10),
        _genome("duplicate", Label.HOST, "TGCA" * 10),
    )
    with pytest.raises(InputError, match="genome IDs must be unique"):
        generate_fragments(genomes, fragment_lengths=(4,), fragments_per_genome=2, seed=1)


@pytest.mark.parametrize(
    ("fragment_lengths", "message"),
    [
        ("31", "iterable of positive integers"),
        (31, "iterable of positive integers"),
        ((), "at least one length"),
        ((True,), "only integers"),
        ((1.5,), "only integers"),
        ((0, -1), "must be positive"),
        ((4, 4), "must not contain duplicates"),
    ],
)
def test_generation_rejects_malformed_fragment_lengths(
    fragment_lengths: object,
    message: str,
) -> None:
    genome = _genome("valid", Label.VIRUS, "ACGT" * 10)
    with pytest.raises(InputError, match=message):
        generate_fragments(
            (genome,),
            fragment_lengths=fragment_lengths,  # type: ignore[arg-type]
            fragments_per_genome=4,
            seed=1,
        )


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"seed": True}, "seed must be an integer"),
        ({"fragments_per_genome": True}, "fragments_per_genome must be an integer"),
        ({"strand_mode": "reverse"}, "strand_mode"),
        ({"max_ambiguous_fraction": True}, "finite number"),
        ({"max_ambiguous_fraction": float("nan")}, "finite number"),
        ({"max_ambiguous_fraction": 1.1}, "finite number"),
    ],
)
def test_generation_option_types_and_bounds_are_checked(
    updates: dict[str, object],
    message: str,
) -> None:
    genome = _genome("valid", Label.VIRUS, "ACGT" * 10)
    arguments: dict[str, object] = {
        "fragment_lengths": (4,),
        "fragments_per_genome": 2,
        "seed": 1,
    }
    arguments.update(updates)
    with pytest.raises(InputError, match=message):
        generate_fragments((genome,), **arguments)  # type: ignore[arg-type]


def test_generation_is_independent_of_genome_contig_and_length_order() -> None:
    first = _genome("genome-a", Label.VIRUS, "ACGT" * 20, "TGCATGCA" * 8)
    second = _genome("genome-b", Label.HOST, "GATTACA" * 12, "CCGGAATT" * 9)
    reordered_first = Genome(
        genome_id=first.genome_id,
        label=first.label,
        contigs=tuple(reversed(first.contigs)),
    )
    reordered_second = Genome(
        genome_id=second.genome_id,
        label=second.label,
        contigs=tuple(reversed(second.contigs)),
    )

    expected = generate_fragments(
        (first, second),
        fragment_lengths=(7, 11),
        fragments_per_genome=17,
        seed=9182,
    )
    reordered = generate_fragments(
        (reordered_second, reordered_first),
        fragment_lengths=(11, 7),
        fragments_per_genome=17,
        seed=9182,
    )

    assert reordered == expected
    assert [fragment.ordinal for fragment in expected] == list(range(len(expected)))
    assert (
        generate_fragments(
            (first, second),
            fragment_lengths=(7, 11),
            fragments_per_genome=17,
            seed=9183,
        )
        != expected
    )


def test_fragment_lengths_are_balanced_within_every_genome() -> None:
    genomes = (
        _genome("v-one", Label.VIRUS, "ACGT" * 100),
        _genome("h-one", Label.HOST, "TGCA" * 100),
    )
    fragments = generate_fragments(
        genomes,
        fragment_lengths=(5, 10, 15),
        fragments_per_genome=10,
        seed=7,
        strand_mode="forward",
    )

    counts: dict[str, Counter[int]] = defaultdict(Counter)
    for fragment in fragments:
        counts[fragment.genome_id][fragment.length] += 1

    assert len(fragments) == 20
    for genome in genomes:
        genome_counts = counts[genome.genome_id]
        assert set(genome_counts) == {5, 10, 15}
        assert sum(genome_counts.values()) == 10
        assert max(genome_counts.values()) - min(genome_counts.values()) <= 1


def test_every_genome_must_support_every_length_before_sampling() -> None:
    short = _genome("short-genome", Label.VIRUS, "ACGT")
    long = _genome("long-genome", Label.HOST, "ACGT" * 20)

    with pytest.raises(InputError) as captured:
        generate_fragments(
            (long, short),
            fragment_lengths=(3, 5),
            fragments_per_genome=4,
            seed=1,
        )

    message = str(captured.value)
    assert "short-genome" in message
    assert "length 5" in message
    assert "within one contig" in message


def test_every_requested_length_receives_at_least_two_fragments() -> None:
    genome = _genome("genome-one", Label.VIRUS, "ACGT" * 20)

    with pytest.raises(InputError, match="at least 6"):
        generate_fragments(
            (genome,),
            fragment_lengths=(4, 5, 6),
            fragments_per_genome=2,
            seed=1,
        )


def test_coordinates_are_uniformly_weighted_across_eligible_contigs() -> None:
    genome = Genome(
        genome_id="weighted-genome",
        label=Label.HOST,
        contigs=(
            Contig(sequence_id="short-contig", sequence="A" * 10),
            Contig(sequence_id="long-contig", sequence="C" * 30),
        ),
    )
    fragments = generate_fragments(
        (genome,),
        fragment_lengths=(5,),
        fragments_per_genome=6_400,
        seed=781,
        strand_mode="forward",
    )

    observed = Counter(fragment.sequence_id for fragment in fragments)
    # Six starts are available on the short contig and 26 on the long one.
    expected_short_fraction = 6 / 32
    observed_short_fraction = observed["short-contig"] / len(fragments)
    assert observed_short_fraction == pytest.approx(expected_short_fraction, abs=0.025)


def test_sampling_is_with_replacement_and_never_crosses_contigs() -> None:
    genome = Genome(
        genome_id="two-contig-genome",
        label=Label.VIRUS,
        contigs=(
            Contig(sequence_id="a-contig", sequence="AAAAAA"),
            Contig(sequence_id="c-contig", sequence="CCCCCC"),
        ),
    )
    fragments = generate_fragments(
        (genome,),
        fragment_lengths=(5,),
        fragments_per_genome=100,
        seed=81,
        strand_mode="forward",
    )

    coordinates = [(fragment.sequence_id, fragment.start, fragment.end) for fragment in fragments]
    assert len(set(coordinates)) < len(coordinates)
    assert {fragment.sequence for fragment in fragments} == {"AAAAA", "CCCCC"}
    assert all(fragment.end - fragment.start == 5 for fragment in fragments)


def test_circular_fragments_wrap_origin_with_unwrapped_truth_coordinates() -> None:
    genome = Genome(
        genome_id="circular-genome",
        label=Label.VIRUS,
        contigs=(Contig("circular-contig", "AACGT", topology="circular"),),
    )
    fragments = generate_fragments(
        (genome,),
        fragment_lengths=(4,),
        fragments_per_genome=40,
        seed=44,
        strand_mode="forward",
    )

    wrapped = [fragment for fragment in fragments if fragment.end > genome.contigs[0].length]
    assert wrapped
    for fragment in wrapped:
        expected = (genome.contigs[0].sequence * 2)[fragment.start : fragment.end]
        assert fragment.sequence == expected


def test_generated_cross_label_exact_content_is_rejected() -> None:
    genomes = (
        _genome("a-virus", Label.VIRUS, "ACGTTGCA"),
        _genome("z-host", Label.HOST, "ACGTTGCA"),
    )

    with pytest.raises(InputError, match="contradictory virus/host labels"):
        generate_fragments(
            genomes,
            fragment_lengths=(8,),
            fragments_per_genome=2,
            seed=9,
            strand_mode="forward",
        )


def test_fragment_identifier_collision_stops_generation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    genome = _genome("collision-genome", Label.VIRUS, "ACGT" * 10)
    monkeypatch.setattr(fragment_module, "_opaque_fragment_id", lambda **_kwargs: "frag-fixed")

    with pytest.raises(IntegrityError, match="fragment-ID collision"):
        generate_fragments(
            (genome,),
            fragment_lengths=(4,),
            fragments_per_genome=2,
            seed=3,
        )


def test_both_strands_are_default_and_truth_stays_on_forward_source() -> None:
    genome = _genome("strand-genome", Label.VIRUS, "AAAACCCCGGGGTTTT" * 10)
    fragments = generate_fragments(
        (genome,),
        fragment_lengths=(13,),
        fragments_per_genome=200,
        seed=121,
    )

    assert {fragment.strand for fragment in fragments} == {"+", "-"}
    for fragment in fragments:
        source = _source_sequence(fragment, (genome,))
        expected = source if fragment.strand == "+" else reverse_complement(source)
        assert fragment.sequence == expected
        assert fragment.start >= 0
        assert fragment.end <= genome.contigs[0].length


def test_forward_mode_preserves_source_sequence() -> None:
    genome = _genome("forward-genome", Label.HOST, "ACGTTGCA" * 20)
    fragments = generate_fragments(
        (genome,),
        fragment_lengths=(9,),
        fragments_per_genome=50,
        seed=53,
        strand_mode="forward",
    )

    assert {fragment.strand for fragment in fragments} == {"+"}
    assert all(fragment.sequence == _source_sequence(fragment, (genome,)) for fragment in fragments)


def test_ambiguity_rejection_never_returns_a_silent_shortfall() -> None:
    genome = _genome("ambiguous-genome", Label.VIRUS, "N" * 50)

    with pytest.raises(InputError) as captured:
        generate_fragments(
            (genome,),
            fragment_lengths=(10,),
            fragments_per_genome=3,
            seed=19,
            max_ambiguous_fraction=0.0,
        )

    message = str(captured.value)
    assert "fragment 1/3" in message
    assert "10,000 random-coordinate attempts" in message
    assert "raise max_ambiguous_fraction" in message


def test_ambiguity_limit_is_inclusive() -> None:
    genome = _genome("threshold-genome", Label.HOST, "AAAANAAAAA")
    fragments = generate_fragments(
        (genome,),
        fragment_lengths=(10,),
        fragments_per_genome=2,
        seed=3,
        strand_mode="forward",
        max_ambiguous_fraction=0.1,
    )

    assert len(fragments) == 2
    assert all(fragment.sequence == "AAAANAAAAA" for fragment in fragments)


def test_fragment_ids_are_unique_opaque_and_label_free() -> None:
    virus = _genome("virus-secret-source", Label.VIRUS, "ACGT" * 40)
    host = _genome("host-secret-source", Label.HOST, "TGCA" * 40)
    fragments = generate_fragments(
        (virus, host),
        fragment_lengths=(20,),
        fragments_per_genome=40,
        seed=991,
    )

    ids = [fragment.fragment_id for fragment in fragments]
    assert len(ids) == len(set(ids))
    assert all(fragment_id.startswith("frag-") and len(fragment_id) == 37 for fragment_id in ids)
    assert all("virus" not in fragment_id and "host" not in fragment_id for fragment_id in ids)
    assert all("secret" not in fragment_id and "source" not in fragment_id for fragment_id in ids)


def test_random_split_is_per_genome_disjoint_and_input_order_independent() -> None:
    genomes = (
        _genome("g-a-virus", Label.VIRUS, "ACGT" * 50),
        _genome("g-b-host", Label.HOST, "TGCA" * 50),
        _genome("g-c-virus", Label.VIRUS, "GATTACA" * 30),
        _genome("g-d-host", Label.HOST, "CCGGAATT" * 30),
    )
    fragments = generate_fragments(
        genomes,
        fragment_lengths=(15,),
        fragments_per_genome=12,
        seed=808,
    )

    train, test = split_fragments_random(fragments, test_fraction=0.25, seed=91)
    reversed_train, reversed_test = split_fragments_random(
        reversed(fragments), test_fraction=0.25, seed=91
    )

    assert (reversed_train, reversed_test) == (train, test)
    train_ids = {fragment.fragment_id for fragment in train}
    test_ids = {fragment.fragment_id for fragment in test}
    assert train_ids.isdisjoint(test_ids)
    assert train_ids | test_ids == {fragment.fragment_id for fragment in fragments}
    for genome in genomes:
        assert sum(fragment.genome_id == genome.genome_id for fragment in train) == 9
        assert sum(fragment.genome_id == genome.genome_id for fragment in test) == 3

    # Final semantic shuffles prevent class/source blocks in serialized output.
    assert [fragment.genome_id for fragment in train] != sorted(
        fragment.genome_id for fragment in train
    )
    assert [fragment.label for fragment in train] != sorted(
        (fragment.label for fragment in train), key=str
    )


def test_random_split_clamps_each_genome_to_nonempty_partitions() -> None:
    genomes = (
        _genome("two-virus", Label.VIRUS, "ACGT" * 20),
        _genome("two-host", Label.HOST, "TGCA" * 20),
    )
    fragments = generate_fragments(
        genomes,
        fragment_lengths=(8,),
        fragments_per_genome=2,
        seed=2,
    )

    train, test = split_fragments_random(fragments, test_fraction=0.99, seed=3)

    assert Counter(fragment.genome_id for fragment in train) == Counter(
        {"two-virus": 1, "two-host": 1}
    )
    assert Counter(fragment.genome_id for fragment in test) == Counter(
        {"two-virus": 1, "two-host": 1}
    )


def test_random_split_rejects_singletons_and_duplicate_ids() -> None:
    genome = _genome("singleton", Label.VIRUS, "ACGT" * 10)
    singleton = Fragment(
        "frag-singleton",
        genome.contigs[0].sequence[:6],
        genome.label,
        genome.genome_id,
        genome.contigs[0].sequence_id,
        0,
        6,
        "+",
        0,
    )

    with pytest.raises(InputError, match="at least two fragments"):
        split_fragments_random((singleton,), test_fraction=0.2, seed=5)
    with pytest.raises(InputError, match="fragment IDs must be unique"):
        split_fragments_random((singleton, singleton), test_fraction=0.2, seed=5)


@pytest.mark.parametrize("test_fraction", [0.0, 1.0, -0.1, float("nan"), float("inf")])
def test_random_split_rejects_invalid_fraction(test_fraction: float) -> None:
    genome = _genome("fraction-genome", Label.HOST, "ACGT" * 10)
    fragments = generate_fragments(
        (genome,),
        fragment_lengths=(6,),
        fragments_per_genome=2,
        seed=6,
    )

    with pytest.raises(InputError, match="greater than 0 and less than 1"):
        split_fragments_random(fragments, test_fraction=test_fraction, seed=7)


@pytest.mark.parametrize(
    ("fragments", "message"),
    [
        ("not-fragments", "iterable of Fragment"),
        (42, "iterable of Fragment"),
        ((), "at least one fragment"),
        ((object(),), "only Fragment"),
    ],
)
def test_random_split_rejects_malformed_fragment_collections(
    fragments: object,
    message: str,
) -> None:
    with pytest.raises(InputError, match=message):
        split_fragments_random(
            fragments,  # type: ignore[arg-type]
            test_fraction=0.2,
            seed=7,
        )


def test_random_split_rejects_seed_and_fraction_types() -> None:
    genome = _genome("typed-split", Label.HOST, "ACGT" * 10)
    fragments = generate_fragments(
        (genome,),
        fragment_lengths=(4,),
        fragments_per_genome=2,
        seed=6,
    )

    with pytest.raises(InputError, match="seed must be an integer"):
        split_fragments_random(fragments, test_fraction=0.2, seed=True)
    with pytest.raises(InputError, match="test_fraction must be a finite number"):
        split_fragments_random(fragments, test_fraction=True, seed=7)
