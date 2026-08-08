"""Versioned JSON Schema resource, contract, and CLI tests."""

from __future__ import annotations

import json
import random
from datetime import date
from pathlib import Path

import pytest
from jsonschema.exceptions import SchemaError

import chimera.schema_resources as schema_module
from chimera.cli import main
from chimera.config import BenchmarkConfig, SplitKind
from chimera.errors import ConfigurationError
from chimera.models import Contig, Fragment, Genome, GenomeMetadata, Label
from chimera.output import (
    REFERENCE_COLUMNS,
    SEQUENCE_COLUMNS,
    TRUTH_COLUMNS,
    reference_rows,
    sequence_rows,
    truth_rows,
)
from chimera.pipeline import (
    ASSIGNMENT_COLUMNS,
    EXCLUSION_COLUMNS,
    _assignment_row,
    generate_benchmark,
)
from chimera.schema_resources import (
    JSON_SCHEMA_NAMES,
    SchemaName,
    SchemaValidationError,
    load_schema,
    schema_filename,
    validate_instance,
    validate_packaged_schemas,
)
from chimera.splits import SplitAssignment, SplitPartition


def _partition_stats() -> dict[str, object]:
    return {
        "records": 2,
        "bases": 62,
        "gc_fraction": 0.5,
        "ambiguous_fraction": 0.0,
        "records_by_label": {"host": 1, "virus": 1},
        "records_by_label_and_length": {"host:31": 1, "virus:31": 1},
        "records_by_length": {"31": 2},
        "source_genomes": 2,
        "records_by_genome": {"host-1": 1, "virus-1": 1},
    }


def _dna(seed: int, *, length: int = 240) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def _split_manifest() -> dict[str, object]:
    return {
        "schema": "urn:chimera:split-manifest:2",
        "protocol": "genome",
        "protocol_id": "genome",
        "parameters": {"test_fraction": 0.5},
        "validation": {
            "status": "pass",
            "diagnostic_only": False,
            "fragment_id_overlap": 0,
            "exact_fragment_content_overlap": 0,
            "test_fragments_with_coordinate_overlap": 0,
            "source_genome_overlap": 0,
            "source_content_hash_overlap": 0,
            "train_source_genomes": 2,
            "test_source_genomes": 2,
        },
        "train": _partition_stats(),
        "test": _partition_stats(),
        "truth_rows": {"train": 2, "test": 2},
        "excluded_genomes": 0,
    }


def _bundle_manifest() -> dict[str, object]:
    virus_digest = "a" * 64
    host_digest = "b" * 64
    return {
        "schema": "urn:chimera:benchmark-bundle:2",
        "tool": {
            "name": "CHIMERA",
            "version": "2.0.0",
            "software_content_sha256": "c" * 64,
            "git_revision": "unknown",
            "git_dirty": None,
        },
        "data_model": {
            "alphabet": "IUPAC DNA",
            "coordinate_system": "0-based-half-open",
            "coordinate_systems": {
                "linear": "0-based-half-open",
                "circular": "0-based-half-open-circular",
            },
            "coordinate_semantics": {
                "linear": "source_start is inclusive and source_end is exclusive",
                "circular": "coordinates are an unwrapped forward-source interval",
            },
            "fragment_headers": "opaque label-free identifiers",
            "grouping": "canonical_topology_aware_genome_sha256_v2",
            "synthetic": True,
        },
        "randomness": {
            "master_seed": 42,
            "algorithm": "Python random.Random",
            "seed_derivation": "chimera.seed.v1",
            "python_implementation": "CPython",
            "python_version": "3.12.0",
        },
        "references": {
            "count": 4,
            "inputs": [
                {
                    "content_id": f"sha256:{virus_digest}",
                    "role": "reference_fasta",
                    "sha256": virus_digest,
                },
                {
                    "content_id": f"sha256:{host_digest}",
                    "role": "reference_fasta",
                    "sha256": host_digest,
                },
            ],
            "preflight_exclusions": 0,
        },
        "splits": {"genome": _split_manifest()},
    }


def test_all_packaged_schemas_are_valid_versioned_resources() -> None:
    assert tuple(name.value for name in SchemaName) == JSON_SCHEMA_NAMES
    validate_packaged_schemas()

    schema_versions = {SchemaName.BUNDLE: 2, SchemaName.SPLIT: 2}
    for name in SchemaName:
        document = load_schema(name)
        assert document["$schema"] == "https://json-schema.org/draft/2020-12/schema"
        assert str(document["$id"]).startswith("urn:chimera:schema:")
        expected_version = schema_versions.get(name, 1)
        assert str(document["$id"]).endswith(f":{expected_version}")


def test_schema_loads_are_fresh_and_unknown_names_are_actionable() -> None:
    first = load_schema("truth-row")
    first["title"] = "mutated by caller"
    assert load_schema("truth-row")["title"] != "mutated by caller"

    with pytest.raises(ConfigurationError, match="Unknown JSON Schema"):
        load_schema("not-a-schema")

    assert schema_filename(SchemaName.BUNDLE) == "benchmark-bundle.schema.json"
    assert schema_filename("truth-row") == "truth-row.schema.json"
    with pytest.raises(ConfigurationError, match="Unknown JSON Schema"):
        schema_filename("not-a-schema")


def test_row_schema_required_fields_track_live_serializer_constants() -> None:
    contracts = {
        "truth-row": TRUTH_COLUMNS,
        "reference-row": REFERENCE_COLUMNS,
        "sequence-row": SEQUENCE_COLUMNS,
        "assignment-row": ASSIGNMENT_COLUMNS,
        "exclusion-row": EXCLUSION_COLUMNS,
    }
    for name, columns in contracts.items():
        required = load_schema(name)["required"]
        assert required == list(columns)


def test_live_reference_sequence_truth_and_assignment_rows_validate(tmp_path: Path) -> None:
    contig = Contig(
        "NC_000001.1",
        "ACGT" * 20,
        source_path=tmp_path / "virus.fna",
        accession_version="NC_000001.1",
        release_date=date(2020, 1, 2),
        metadata_extra=(("isolate", "example"),),
    )
    genome = Genome(
        "virus-1",
        Label.VIRUS,
        (contig,),
        GenomeMetadata(
            release_date=date(2020, 1, 2),
            accession_version="NC_000001.1",
            taxonomy=(("family", "Example;viridae"),),
            extra=(("collection", "synthetic"),),
        ),
    )
    fragment = Fragment(
        "frag-0123456789abcdef0123456789abcdef",
        "ACGT",
        Label.VIRUS,
        "virus-1",
        "NC_000001.1",
        0,
        4,
        "+",
        0,
    )
    truth = truth_rows((fragment,), partition="test", genomes={genome.genome_id: genome})[0]
    input_digest = "e" * 64
    source_digests = {contig.source_path.resolve(): input_digest}

    validate_instance(
        reference_rows((genome,), source_digests=source_digests)[0],
        "reference-row",
    )
    validate_instance(
        sequence_rows((genome,), source_digests=source_digests)[0],
        "sequence-row",
    )
    validate_instance(truth, "truth-row")

    assignment = SplitAssignment(
        genome_id=genome.genome_id,
        label=genome.label,
        partition=SplitPartition.TEST,
        reason="label_stratified_genome_holdout",
        group_id=f"sha256:{genome.digest}",
    )
    validate_instance(_assignment_row(assignment), "assignment-row")


def test_exclusion_metadata_and_resolved_config_rows_validate(tmp_path: Path) -> None:
    digest = "b" * 64
    exclusion = dict.fromkeys(EXCLUSION_COLUMNS, "")
    exclusion.update(
        {
            "genome_id": "virus-duplicate",
            "label": "virus",
            "split": "reference_preflight",
            "reason": "same_class_content_duplicate",
            "duplicate_of": "virus-primary",
            "source_sha256": digest,
        }
    )
    validate_instance(exclusion, "exclusion-row")
    validate_instance(
        {
            "sequence_id": "NC_000001.1",
            "genome_id": "virus-1",
            "label": "virus",
            "accession_version": "",
            "release_date": "2020-01-02",
            "topology": "circular",
            "family": "Exampleviridae",
            "custom_provenance": "value",
        },
        "metadata-row",
    )

    config = BenchmarkConfig(
        virus_paths=(tmp_path / "virus.fna",),
        host_paths=(tmp_path / "host.fna",),
        output_dir=tmp_path / "bundle",
        splits=(SplitKind.GENOME,),
        fragment_lengths=(31,),
        fragments_per_genome=2,
    )
    resolved = config.as_manifest_dict()
    resolved.pop("overwrite")
    resolved["virus_paths"] = [f"sha256:{'a' * 64}"]
    resolved["host_paths"] = [f"sha256:{'b' * 64}"]
    resolved["output_dir"] = "bundle"
    resolved = json.loads(json.dumps(resolved))
    validate_instance(resolved, "resolved-config")


@pytest.mark.parametrize(
    "value",
    ["0", "0.0", "0.95", "0.12345678901234568", "1e-07", "1.0", "1"],
)
def test_fractional_tsv_evidence_accepts_shortest_round_trip_strings(value: str) -> None:
    row = dict.fromkeys(ASSIGNMENT_COLUMNS, "")
    row.update(
        {
            "genome_id": "virus-test",
            "group_id": f"sha256:{'c' * 64}",
            "label": "virus",
            "partition": "test",
            "reason": "similarity_filtered_test",
            "max_train_similarity": value,
        }
    )
    validate_instance(row, "assignment-row")


@pytest.mark.parametrize("value", ["-0.1", "1.0000000001", "NaN", "inf"])
def test_fractional_tsv_evidence_rejects_values_outside_unit_interval(value: str) -> None:
    row = dict.fromkeys(ASSIGNMENT_COLUMNS, "")
    row.update(
        {
            "genome_id": "virus-test",
            "group_id": f"sha256:{'d' * 64}",
            "label": "virus",
            "partition": "test",
            "reason": "similarity_filtered_test",
            "max_train_similarity": value,
        }
    )
    with pytest.raises(SchemaValidationError, match="max_train_similarity"):
        validate_instance(row, "assignment-row")


def test_bundle_validation_resolves_packaged_split_schema_offline() -> None:
    validate_instance(_bundle_manifest(), "bundle")

    corrupted = _bundle_manifest()
    corrupted["splits"] = {"genome": {}}
    with pytest.raises(SchemaValidationError, match="bundle schema violation"):
        validate_instance(corrupted, "bundle")


def test_generated_json_documents_satisfy_the_published_contracts(tmp_path: Path) -> None:
    virus = tmp_path / "virus.fna"
    host = tmp_path / "host.fna"
    metadata = tmp_path / "metadata.tsv"
    virus.write_text(f">v-old\n{_dna(1)}\n>v-new\n{_dna(2)}\n", encoding="utf-8")
    host.write_text(f">h-old\n{_dna(3)}\n>h-new\n{_dna(4)}\n", encoding="utf-8")
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\trelease_date\tfamily\n"
        "v-old\tv-old\tvirus\t2018-01-01\tAlpha\n"
        "v-new\tv-new\tvirus\t2022-01-01\tBeta\n"
        "h-old\th-old\thost\t2017-01-01\t\n"
        "h-new\th-new\thost\t2023-01-01\t\n",
        encoding="utf-8",
    )
    config = BenchmarkConfig(
        virus_paths=(virus,),
        host_paths=(host,),
        metadata_path=metadata,
        output_dir=tmp_path / "bundle",
        fragment_lengths=(31,),
        fragments_per_genome=4,
        similarity_k=15,
        sketch_size=100,
        seed=1729,
    )
    result = generate_benchmark(config)

    manifest = json.loads((result.output_dir / "manifest.json").read_text(encoding="utf-8"))
    resolved = json.loads((result.output_dir / "resolved-config.json").read_text(encoding="utf-8"))
    receipts = manifest["references"]["inputs"]
    assert all(receipt["content_id"] == f"sha256:{receipt['sha256']}" for receipt in receipts)
    assert resolved["output_dir"] == "bundle"
    assert all(value.startswith("sha256:") for value in resolved["virus_paths"])
    assert all(value.startswith("sha256:") for value in resolved["host_paths"])
    assert str(tmp_path) not in json.dumps({"manifest": manifest, "resolved": resolved})
    validate_instance(manifest, "bundle")
    validate_instance(resolved, "resolved-config")
    for kind in SplitKind:
        split_manifest = json.loads(
            (result.output_dir / kind.directory_name / "split.json").read_text(encoding="utf-8")
        )
        validate_instance(split_manifest, "split")


def test_schema_cli_preserves_headers_and_prints_every_json_schema(capsys) -> None:
    header_contracts = {
        "truth": TRUTH_COLUMNS,
        "references": REFERENCE_COLUMNS,
    }
    for name, expected in header_contracts.items():
        assert main(["schema", name]) == 0
        assert tuple(capsys.readouterr().out.strip().split("\t")) == expected

    for name in JSON_SCHEMA_NAMES:
        assert main(["schema", name]) == 0
        document = json.loads(capsys.readouterr().out)
        assert document == load_schema(name)


def test_schema_resource_read_failures_are_actionable(monkeypatch):
    class BrokenResource:
        def is_file(self):
            raise OSError("unreadable")

        def __str__(self):
            return "broken-resource"

    monkeypatch.setattr(
        schema_module,
        "_schema_candidates",
        lambda _filename: (BrokenResource(),),
    )
    with pytest.raises(ConfigurationError, match=r"unavailable.*unreadable"):
        load_schema("bundle")


def test_malformed_and_nonobject_schema_resources_are_rejected(monkeypatch):
    monkeypatch.setattr(schema_module, "_read_schema_text", lambda _name: "{")
    with pytest.raises(ConfigurationError, match="malformed"):
        load_schema("bundle")

    monkeypatch.setattr(schema_module, "_read_schema_text", lambda _name: "[]")
    with pytest.raises(ConfigurationError, match="not a JSON object"):
        load_schema("bundle")

    assert schema_module._is_json_value({"nested": object()}) is False


def test_offline_registry_requires_every_schema_identifier(monkeypatch):
    monkeypatch.setattr(
        schema_module,
        "load_schema",
        lambda _name: {"$schema": "https://json-schema.org/draft/2020-12/schema"},
    )
    with pytest.raises(ConfigurationError, match="has no nonempty \\$id"):
        schema_module._offline_registry()


def test_validate_instance_wraps_an_invalid_packaged_schema(monkeypatch):
    def invalid(_document):
        raise SchemaError("invalid schema")

    monkeypatch.setattr(schema_module.Draft202012Validator, "check_schema", invalid)
    with pytest.raises(ConfigurationError, match=r"Packaged JSON Schema.*invalid"):
        validate_instance({}, "bundle")


def test_validate_packaged_schemas_wraps_schema_errors(monkeypatch):
    real_check = schema_module.Draft202012Validator.check_schema
    calls = 0

    def fail_once(document):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise SchemaError("invalid packaged contract")
        return real_check(document)

    monkeypatch.setattr(schema_module.Draft202012Validator, "check_schema", fail_once)
    with pytest.raises(ConfigurationError, match=r"Packaged JSON Schema.*invalid"):
        validate_packaged_schemas()
