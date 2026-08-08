from __future__ import annotations

from datetime import date
from pathlib import Path

import pytest

from chimera.config import (
    BenchmarkConfig,
    ConfigurationError,
    DuplicatePolicy,
    MissingMetadataPolicy,
    SimilarityBands,
    SplitKind,
    parse_splits,
)


def minimal_mapping() -> dict[str, object]:
    return {
        "virus_paths": ["virus.fna"],
        "host_paths": ["host.fna"],
        "output_dir": "result",
    }


def test_defaults_select_all_protocols_and_resolve_relative_paths(tmp_path):
    config = BenchmarkConfig.from_mapping(minimal_mapping(), base_dir=tmp_path)

    assert config.splits == tuple(SplitKind)
    assert config.virus_paths == ((tmp_path / "virus.fna").resolve(),)
    assert config.output_dir == (tmp_path / "result").resolve()


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("random-fragment", SplitKind.RANDOM),
        ("genome-level", SplitKind.GENOME),
        ("similarity-filtered", SplitKind.SIMILARITY),
        ("temporal", SplitKind.TEMPORAL),
        ("taxonomic-holdout", SplitKind.TAXONOMY),
    ],
)
def test_split_aliases_are_intuitive(value, expected):
    assert SplitKind.parse(value) is expected


def test_split_directory_names_are_stable() -> None:
    assert {kind: kind.directory_name for kind in SplitKind} == {
        SplitKind.RANDOM: "random_fragment",
        SplitKind.GENOME: "genome_holdout",
        SplitKind.SIMILARITY: "similarity_filtered",
        SplitKind.TEMPORAL: "temporal_holdout",
        SplitKind.TAXONOMY: "taxonomic_holdout",
    }


def test_unknown_and_empty_split_selections_are_actionable() -> None:
    with pytest.raises(ConfigurationError, match=r"Unknown split.*choose from"):
        SplitKind.parse("not-a-protocol")
    with pytest.raises(ConfigurationError, match="At least one split"):
        parse_splits(["", "  ,  "])
    with pytest.raises(ConfigurationError, match="splits must be"):
        parse_splits(42)


def test_split_list_is_deduplicated_in_protocol_order():
    assert parse_splits("taxonomy,random-fragment,random,genome") == (
        SplitKind.RANDOM,
        SplitKind.GENOME,
        SplitKind.TAXONOMY,
    )


def test_toml_uses_file_directory_for_paths(tmp_path):
    config_path = tmp_path / "benchmark.toml"
    config_path.write_text(
        """
[benchmark]
schema_version = 1
virus_paths = ["data/virus.fna"]
host_paths = ["data/host.fna"]
output_dir = "output"
splits = ["genome", "temporal"]
temporal_cutoff = 2020-01-01
fragment_lengths = [150, 500]
""".strip(),
        encoding="utf-8",
    )

    config = BenchmarkConfig.from_toml(config_path)

    assert config.temporal_cutoff == date(2020, 1, 1)
    assert config.splits == (SplitKind.GENOME, SplitKind.TEMPORAL)
    assert config.fragment_lengths == (150, 500)
    assert config.virus_paths[0] == (tmp_path / "data/virus.fna").resolve()


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"virus_paths": []}, "virus"),
        ({"test_fraction": 1.0}, "test_fraction"),
        ({"fragment_lengths": [0]}, "fragment_lengths"),
        ({"similarity_k": 20}, "similarity_k"),
        ({"sketch_size": 99}, "sketch_size"),
        ({"unknown": True}, "Unknown configuration"),
    ],
)
def test_invalid_configurations_fail_with_actionable_message(tmp_path, updates, message):
    values = minimal_mapping()
    values.update(updates)
    with pytest.raises((ConfigurationError, ValueError), match=message):
        BenchmarkConfig.from_mapping(values, base_dir=tmp_path)


@pytest.mark.parametrize(
    ("updates", "field_name"),
    [
        ({"similarity_bands": [0.9, 0.7, 0.3]}, "similarity_bands"),
        ({"similarity_bands": {"high": "near"}}, "similarity_bands"),
        ({"fragment_lengths": 31.5}, "fragment_lengths"),
        ({"fragment_lengths": None}, "fragment_lengths"),
        ({"holdout_taxa": 7}, "holdout_taxa"),
        ({"holdout_taxa": None}, "holdout_taxa"),
        ({"output_dir": ""}, "output_dir"),
        ({"output_dir": Path()}, "output_dir"),
        ({"output_dir": []}, "output_dir"),
    ],
)
def test_malformed_field_shapes_raise_configuration_error(tmp_path, updates, field_name):
    values = minimal_mapping()
    values.update(updates)

    with pytest.raises(ConfigurationError, match=field_name):
        BenchmarkConfig.from_mapping(values, base_dir=tmp_path)


def test_similarity_bands_are_mutually_exclusive():
    bands = SimilarityBands(high=0.9, moderate=0.7, low=0.3)
    assert bands.classify(0.95) == "high_similarity"
    assert bands.classify(0.8) == "moderate_similarity"
    assert bands.classify(0.4) == "low_similarity"
    assert bands.classify(0.1) == "distant_detectable"
    assert bands.classify(None) == "no_detectable_match"


@pytest.mark.parametrize(
    "bands",
    [
        {"high": float("nan")},
        {"high": 0.7, "moderate": 0.8, "low": 0.3},
        {"unexpected": 0.5},
    ],
)
def test_similarity_band_validation_rejects_nonfinite_ordered_and_unknown_values(
    tmp_path: Path,
    bands: dict[str, object],
) -> None:
    values = minimal_mapping() | {"similarity_bands": bands}
    with pytest.raises(ConfigurationError, match=r"Similarity band|similarity_bands"):
        BenchmarkConfig.from_mapping(values, base_dir=tmp_path)


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"host_paths": []}, "host"),
        ({"fragment_lengths": [31, 31]}, "duplicates"),
        ({"fragment_lengths": [31, 61], "fragments_per_genome": 3}, "at least 4"),
        ({"fragment_lengths": ["not-an-integer"]}, "fragment_lengths"),
        ({"strand_mode": "reverse"}, "strand_mode"),
        ({"max_ambiguous_fraction": 1.1}, "max_ambiguous_fraction"),
        ({"taxonomy_rank": "viral family"}, "taxonomy_rank"),
        ({"auto_holdout_count": 0}, "auto_holdout_count"),
        ({"max_train_similarity": -0.1}, "max_train_similarity"),
        ({"min_similarity_coverage": 1.1}, "min_similarity_coverage"),
        (
            {
                "max_train_similarity": 0.89,
                "similarity_bands": {"high": 0.9, "moderate": 0.7, "low": 0.3},
            },
            "high similarity-band",
        ),
        ({"schema_version": 2}, "Unsupported config schema_version"),
        ({"temporal_cutoff": "not-a-date"}, "ISO date"),
        ({"virus_paths": ["same.fna", "same.fna"]}, "same path"),
        ({"virus_paths": object()}, "virus_paths"),
        ({"output_dir": object()}, "output_dir"),
    ],
)
def test_additional_configuration_error_paths_are_actionable(
    tmp_path: Path,
    updates: dict[str, object],
    message: str,
) -> None:
    values = minimal_mapping() | updates
    with pytest.raises(ConfigurationError, match=message):
        BenchmarkConfig.from_mapping(values, base_dir=tmp_path)


def test_scalar_and_string_fields_are_normalized(tmp_path: Path) -> None:
    absolute_virus = (tmp_path / "absolute-virus.fna").resolve()
    config = BenchmarkConfig.from_mapping(
        minimal_mapping()
        | {
            "virus_paths": absolute_virus,
            "fragment_lengths": 31,
            "holdout_taxa": " Alpha ",
            "taxonomy_rank": "FAMILY",
        },
        base_dir=tmp_path,
    )

    assert config.virus_paths == (absolute_virus,)
    assert config.fragment_lengths == (31,)
    assert config.holdout_taxa == ("Alpha",)
    assert config.taxonomy_rank == "family"


def test_direct_config_rejects_empty_splits(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="evaluation split"):
        BenchmarkConfig(
            virus_paths=(tmp_path / "virus.fna",),
            host_paths=(tmp_path / "host.fna",),
            output_dir=tmp_path / "bundle",
            splits=(),
        )


def test_toml_read_and_section_shape_errors_are_wrapped(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="Cannot read TOML configuration"):
        BenchmarkConfig.from_toml(tmp_path / "missing.toml")

    malformed = tmp_path / "malformed.toml"
    malformed.write_text("[benchmark\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="Cannot read TOML configuration"):
        BenchmarkConfig.from_toml(malformed)

    non_table = tmp_path / "non-table.toml"
    non_table.write_text('benchmark = "not a table"\n', encoding="utf-8")
    with pytest.raises(ConfigurationError, match=r"\[benchmark\] must be a table"):
        BenchmarkConfig.from_toml(non_table)


def test_manifest_dict_is_json_compatible(tmp_path):
    config = BenchmarkConfig.from_mapping(
        minimal_mapping()
        | {
            "metadata_path": "metadata.tsv",
            "similarity_table": "similarity.tsv",
            "temporal_cutoff": "2020-01-02",
            "duplicate_policy": DuplicatePolicy.DROP,
            "missing_metadata": MissingMetadataPolicy.EXCLUDE,
        },
        base_dir=tmp_path,
    )
    manifest = config.as_manifest_dict()
    assert manifest["schema_version"] == 1
    assert manifest["splits"] == [kind.value for kind in SplitKind]
    assert manifest["duplicate_policy"] == "drop"
    assert manifest["missing_metadata"] == "exclude"
    assert manifest["metadata_path"] == str((tmp_path / "metadata.tsv").resolve())
    assert manifest["similarity_table"] == str((tmp_path / "similarity.tsv").resolve())
    assert manifest["temporal_cutoff"] == "2020-01-02"
