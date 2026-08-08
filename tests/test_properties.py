from __future__ import annotations

import pytest

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given, settings
from hypothesis import strategies as st

from chimera.fragments import generate_fragments
from chimera.models import (
    Contig,
    Genome,
    Label,
    canonical_sequence_hash,
    normalize_iupac_dna,
    reverse_complement,
)
from chimera.splits import genome_holdout

IUPAC = "ACGTRYSWKMBDHVN"


@given(st.text(alphabet=IUPAC, min_size=1, max_size=500))
def test_reverse_complement_is_an_involution(sequence):
    assert reverse_complement(reverse_complement(sequence)) == normalize_iupac_dna(sequence)


@given(st.text(alphabet=IUPAC, min_size=1, max_size=500))
def test_canonical_hash_is_orientation_invariant(sequence):
    assert canonical_sequence_hash(sequence) == canonical_sequence_hash(
        reverse_complement(sequence)
    )


@settings(max_examples=40, deadline=None)
@given(
    sequence=st.text(alphabet="ACGT", min_size=40, max_size=300),
    fragment_length=st.integers(min_value=1, max_value=40),
    seed=st.integers(),
)
def test_generated_fragment_coordinates_reconstruct_the_forward_template(
    sequence, fragment_length, seed
):
    genome = Genome("G", Label.VIRUS, (Contig("G.1", sequence),))
    fragments = generate_fragments(
        [genome],
        fragment_lengths=[fragment_length],
        fragments_per_genome=5,
        seed=seed,
    )
    for fragment in fragments:
        forward = sequence[fragment.start : fragment.end]
        expected = forward if fragment.strand == "+" else reverse_complement(forward)
        assert fragment.sequence == expected
        assert fragment.length == fragment_length


@settings(max_examples=30, deadline=None)
@given(seed=st.integers(), test_fraction=st.floats(min_value=0.05, max_value=0.95))
def test_genome_holdout_remains_disjoint_for_all_seeds_and_fractions(seed, test_fraction):
    genomes = (
        Genome("v1", Label.VIRUS, (Contig("v1.1", "A" * 80 + "C"),)),
        Genome("v2", Label.VIRUS, (Contig("v2.1", "C" * 80 + "G"),)),
        Genome("v3", Label.VIRUS, (Contig("v3.1", "G" * 80 + "T"),)),
        Genome("h1", Label.HOST, (Contig("h1.1", "T" * 80 + "A"),)),
        Genome("h2", Label.HOST, (Contig("h2.1", "AC" * 40 + "G"),)),
        Genome("h3", Label.HOST, (Contig("h3.1", "GT" * 40 + "C"),)),
    )
    plan = genome_holdout(genomes, test_fraction=test_fraction, seed=seed)
    assert plan.train_ids.isdisjoint(plan.test_ids)
    assert {item.label for item in plan.train} == {Label.VIRUS, Label.HOST}
    assert {item.label for item in plan.test} == {Label.VIRUS, Label.HOST}
