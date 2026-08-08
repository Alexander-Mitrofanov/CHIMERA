from __future__ import annotations

import math
from pathlib import Path

import pytest

from chimera.errors import InputError
from chimera.models import Contig, Genome, Label, reverse_complement
from chimera.similarity import (
    GenomeSketch,
    SimilarityHit,
    best_train_matches,
    format_similarity_value,
    mash_identity,
    read_similarity_table,
    similarity_rows,
    sketch_genome,
    sketch_jaccard,
)

TABLE_HEADER = (
    "query_genome_id\treference_genome_id\tsimilarity\tcoverage\tcoverage_definition\tmethod\n"
)


def genome(genome_id: str, sequence: str, label: Label = Label.VIRUS) -> Genome:
    return Genome(genome_id, label, (Contig(f"{genome_id}.1", sequence),))


def test_identical_and_reverse_complement_genomes_have_identity_one():
    left = sketch_genome(genome("a", "ACGT" * 30), k=7, sketch_size=100)
    right = sketch_genome(genome("b", reverse_complement("ACGT" * 30)), k=7, sketch_size=100)

    assert sketch_jaccard(left, right)[0] == 1.0


def test_no_shared_kmers_is_explicitly_undetectable():
    hits = best_train_matches(
        [genome("train", "A" * 100)],
        [genome("test", "C" * 100)],
        k=7,
        sketch_size=100,
    )

    assert hits["test"].similarity is None
    assert hits["test"].reference_genome_id is None


def test_tie_resolution_does_not_depend_on_input_order():
    query = genome("query", "ACGT" * 30)
    train_a = genome("a", "ACGT" * 30)
    train_b = genome("b", "ACGT" * 30)

    first = best_train_matches([train_b, train_a], [query], k=7, sketch_size=100)
    second = best_train_matches([train_a, train_b], [query], k=7, sketch_size=100)

    assert first == second
    assert first["query"].reference_genome_id == "a"


def test_ambiguity_only_genome_cannot_be_sketch_screened():
    with pytest.raises(InputError, match="no unambiguous"):
        sketch_genome(genome("ambiguous", "N" * 100), k=7, sketch_size=100)


@pytest.mark.parametrize("jaccard", [-0.1, 1.1])
def test_mash_identity_rejects_out_of_range_values(jaccard):
    with pytest.raises(ValueError, match="between"):
        mash_identity(jaccard, k=21)


def test_mash_identity_maps_partial_overlap_between_zero_and_one():
    identity = mash_identity(0.5, k=21)

    assert identity is not None
    assert 0.0 < identity < 1.0


def test_explicit_similarity_table_selects_best_train_hit(tmp_path: Path):
    path = tmp_path / "similarity.tsv"
    path.write_text(
        TABLE_HEADER
        + "q\ta\t0.91\t0.88\taligned_fraction_shorter\tskani-0.3.1\n"
        + "q\tb\t0.97\t0.90\taligned_fraction_shorter\tskani-0.3.1\n",
        encoding="utf-8",
    )

    hits = read_similarity_table(path, query_ids={"q"}, reference_ids={"a", "b"})

    assert hits["q"].reference_genome_id == "b"
    assert hits["q"].similarity == 0.97
    assert hits["q"].coverage == 0.90


def test_similarity_table_cannot_hide_qualifying_hit_behind_low_coverage_maximum(tmp_path):
    path = tmp_path / "similarity.tsv"
    path.write_text(
        TABLE_HEADER
        + "q\ta\t0.99\t0.10\taligned_fraction_shorter\tskani-0.3.1\n"
        + "q\tb\t0.98\t0.90\taligned_fraction_shorter\tskani-0.3.1\n",
        encoding="utf-8",
    )

    hit = read_similarity_table(
        path,
        query_ids={"q"},
        reference_ids={"a", "b"},
        max_train_similarity=0.95,
        min_similarity_coverage=0.85,
    )["q"]

    assert hit.reference_genome_id == "a"
    assert hit.similarity == 0.99
    assert hit.coverage == 0.10
    assert hit.strict_gate_reference_genome_id == "b"
    assert hit.strict_gate_similarity == 0.98
    assert hit.strict_gate_coverage == 0.90


def test_similarity_table_requires_explicit_result_for_every_query(tmp_path):
    path = tmp_path / "similarity.tsv"
    path.write_text(
        TABLE_HEADER + "q1\ta\t0.9\t0.9\taligned_fraction_shorter\tskani-0.3.1\n",
        encoding="utf-8",
    )
    with pytest.raises(InputError, match="Cartesian matrix"):
        read_similarity_table(path, query_ids={"q1", "q2"}, reference_ids={"a"})


def test_similarity_table_rejects_unknown_reference(tmp_path):
    path = tmp_path / "similarity.tsv"
    path.write_text(
        TABLE_HEADER + "q\tother\t0.9\t0.9\taligned_fraction_shorter\tskani-0.3.1\n",
        encoding="utf-8",
    )
    with pytest.raises(InputError, match="not a training genome"):
        read_similarity_table(path, query_ids={"q"}, reference_ids={"a"})


def test_similarity_table_ignores_nondetected_pairs_when_collapsing_mixed_rows(tmp_path):
    path = tmp_path / "similarity.tsv"
    path.write_text(
        TABLE_HEADER
        + "q\ta\t\t\taligned_fraction_shorter\tskani-0.3.1\n"
        + "q\tb\t0.72\t0.86\taligned_fraction_shorter\tskani-0.3.1\n",
        encoding="utf-8",
    )

    hit = read_similarity_table(path, query_ids={"q"}, reference_ids={"a", "b"})["q"]

    assert hit.reference_genome_id == "b"
    assert hit.similarity == 0.72
    assert hit.coverage == 0.86
    assert hit.strict_gate_reference_genome_id is None


def test_similarity_table_all_nondetected_pairs_have_no_arbitrary_nearest_id(tmp_path):
    path = tmp_path / "similarity.tsv"
    path.write_text(
        TABLE_HEADER
        + "q\ta\t\t\taligned_fraction_shorter\tskani-0.3.1\n"
        + "q\tb\t\t\taligned_fraction_shorter\tskani-0.3.1\n",
        encoding="utf-8",
    )

    hit = read_similarity_table(path, query_ids={"q"}, reference_ids={"b", "a"})["q"]

    assert hit.reference_genome_id is None
    assert hit.similarity is None
    assert hit.coverage is None
    assert hit.strict_gate_reference_genome_id is None
    assert hit.strict_gate_similarity is None


@pytest.mark.parametrize(("similarity", "coverage"), [("0.8", ""), ("", "0.8")])
def test_similarity_table_rejects_half_blank_detection_evidence(tmp_path, similarity, coverage):
    path = tmp_path / "similarity.tsv"
    path.write_text(
        TABLE_HEADER + f"q\ta\t{similarity}\t{coverage}\taligned_fraction_shorter\tskani-0.3.1\n",
        encoding="utf-8",
    )
    with pytest.raises(InputError, match="both be blank"):
        read_similarity_table(path, query_ids={"q"}, reference_ids={"a"})


def test_similarity_table_rejects_duplicate_pair_rows(tmp_path):
    path = tmp_path / "similarity.tsv"
    path.write_text(
        TABLE_HEADER
        + "q\ta\t\t\taligned_fraction_shorter\tskani-0.3.1\n"
        + "q\ta\t0.8\t0.9\taligned_fraction_shorter\tskani-0.3.1\n",
        encoding="utf-8",
    )

    with pytest.raises(InputError, match="duplicate query/reference pair"):
        read_similarity_table(path, query_ids={"q"}, reference_ids={"a"})


def test_similarity_evidence_serialization_is_lossless():
    similarity = math.nextafter(0.95, 1.0)
    coverage = math.nextafter(0.85, 0.0)
    hit = SimilarityHit(
        query_genome_id="q",
        reference_genome_id="a",
        similarity=similarity,
        coverage=coverage,
        method="skani-0.3.1",
    )

    row = similarity_rows({"q": hit})[0]

    assert float(format_similarity_value(similarity)) == similarity
    assert float(row["similarity"]) == similarity
    assert float(row["coverage"]) == coverage
    assert row["similarity"] != f"{similarity:.8f}"
    assert format_similarity_value(None) == ""


@pytest.mark.parametrize("value", [True, "0.5", object()])
def test_similarity_formatter_rejects_non_numeric_values(value):
    with pytest.raises(TypeError, match="number or None"):
        format_similarity_value(value)


@pytest.mark.parametrize("value", [-0.1, 1.1, math.inf, -math.inf, math.nan])
def test_similarity_formatter_rejects_non_finite_or_out_of_range_values(value):
    with pytest.raises(ValueError, match=r"within \[0, 1\]"):
        format_similarity_value(value)


@pytest.mark.parametrize(("k", "sketch_size"), [(0, 10), (3, 0)])
def test_sketch_rejects_nonpositive_parameters(k, sketch_size):
    with pytest.raises(ValueError, match="positive"):
        sketch_genome(genome("query", "ACGT" * 10), k=k, sketch_size=sketch_size)


def test_jaccard_rejects_mismatched_k_and_handles_empty_approximate_sample():
    left = GenomeSketch("left", frozenset(), 7, 1, 2, "a")
    right = GenomeSketch("right", frozenset(), 9, 1, 2, "b")
    with pytest.raises(ValueError, match="different k"):
        sketch_jaccard(left, right)

    right = GenomeSketch("right", frozenset(), 7, 1, 2, "b")
    assert sketch_jaccard(left, right) == (0.0, 0)


def test_best_train_matches_requires_training_genomes():
    with pytest.raises(ValueError, match="training genome"):
        best_train_matches([], [genome("query", "ACGT" * 10)])


@pytest.mark.parametrize(
    ("keyword", "value"),
    [
        ("max_train_similarity", -0.1),
        ("max_train_similarity", math.nan),
        ("min_similarity_coverage", 1.1),
        ("min_similarity_coverage", math.inf),
    ],
)
def test_similarity_table_rejects_invalid_gate_thresholds(tmp_path, keyword, value):
    options = {keyword: value}
    with pytest.raises(ValueError, match=keyword):
        read_similarity_table(
            tmp_path / "unused.tsv",
            query_ids=["q"],
            reference_ids=["r"],
            **options,
        )


@pytest.mark.parametrize(
    ("queries", "references", "message"),
    [
        (["q", "q"], ["r"], "must not contain duplicates"),
        (["q"], ["r", "r"], "must not contain duplicates"),
        ([], ["r"], "both be nonempty"),
        (["q"], [], "both be nonempty"),
    ],
)
def test_similarity_table_rejects_invalid_expected_id_sets(tmp_path, queries, references, message):
    with pytest.raises(ValueError, match=message):
        read_similarity_table(
            tmp_path / "unused.tsv",
            query_ids=queries,
            reference_ids=references,
        )


def test_similarity_table_wraps_open_failure(tmp_path):
    with pytest.raises(InputError, match="Cannot read similarity table"):
        read_similarity_table(tmp_path / "missing.tsv", query_ids=["q"], reference_ids=["r"])


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ("query_genome_id\nq\n", "missing column"),
        (
            TABLE_HEADER + "other\tr\t0.5\t0.5\taligned_fraction_shorter\tskani-0.3.1\n",
            "unknown/non-test query",
        ),
        (
            TABLE_HEADER + "q\t\t0.5\t0.5\taligned_fraction_shorter\tskani-0.3.1\n",
            "reference_genome_id is required",
        ),
        (
            TABLE_HEADER + "q\tr\tbad\t0.5\taligned_fraction_shorter\tskani-0.3.1\n",
            "similarity must be a number",
        ),
        (
            TABLE_HEADER + "q\tr\t1.1\t0.5\taligned_fraction_shorter\tskani-0.3.1\n",
            "similarity .* outside",
        ),
        (
            TABLE_HEADER + "q\tr\t0.5\tbad\taligned_fraction_shorter\tskani-0.3.1\n",
            "coverage must be a number",
        ),
        (
            TABLE_HEADER + "q\tr\t0.5\t-0.1\taligned_fraction_shorter\tskani-0.3.1\n",
            "coverage .* outside",
        ),
        (
            TABLE_HEADER + "q\tr\t0.5\t0.5\tquery_fraction\tskani-0.3.1\n",
            "coverage_definition",
        ),
        (
            TABLE_HEADER + "q\tr\t0.5\t0.5\taligned_fraction_shorter\tskani\n",
            "tool/version",
        ),
    ],
)
def test_similarity_table_rejects_malformed_rows(tmp_path, payload, message):
    path = tmp_path / "similarity.tsv"
    path.write_text(payload, encoding="utf-8")
    with pytest.raises(InputError, match=message):
        read_similarity_table(path, query_ids=["q"], reference_ids=["r"])
