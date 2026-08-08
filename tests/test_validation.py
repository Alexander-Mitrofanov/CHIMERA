"""Corruption and leakage traps for independent bundle validation."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import random
import shutil
from collections import Counter
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

import pytest

import chimera.validation as validation_module
from chimera.config import BenchmarkConfig, SplitKind
from chimera.errors import IntegrityError
from chimera.fasta import read_fasta, write_fasta
from chimera.models import (
    Contig,
    Topology,
    canonical_sequence_hash,
    deterministic_genome_hash,
)
from chimera.output import write_checksums, write_json, write_tsv
from chimera.pipeline import generate_benchmark
from chimera.schema_resources import SchemaValidationError
from chimera.validation import ValidationReport, validate_bundle


def _dna(seed: int, length: int = 240) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


@pytest.fixture(scope="module")
def bundle_template(tmp_path_factory: pytest.TempPathFactory) -> Path:
    template_root = tmp_path_factory.mktemp("validation-template")
    virus = template_root / "viruses.fna"
    host = template_root / "hosts.fna"
    metadata = template_root / "metadata.tsv"
    virus.write_text(f">v_old\n{_dna(1)}\n>v_new\n{_dna(2)}\n", encoding="utf-8")
    host.write_text(f">h_old\n{_dna(3)}\n>h_new\n{_dna(4)}\n", encoding="utf-8")
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\taccession_version\trelease_date\tfamily\n"
        "v_old\tv_old\tvirus\tNC_000001.1\t2018-01-01\tAlpha\n"
        "v_new\tv_new\tvirus\tNC_000002.1\t2022-01-01\tBeta\n"
        "h_old\th_old\thost\tGCF_000001.1\t2017-01-01\t\n"
        "h_new\th_new\thost\tGCF_000002.1\t2023-01-01\t\n",
        encoding="utf-8",
    )
    output = template_root / "bundle"
    generate_benchmark(
        BenchmarkConfig(
            virus_paths=(virus,),
            host_paths=(host,),
            metadata_path=metadata,
            output_dir=output,
            fragment_lengths=(31,),
            fragments_per_genome=4,
            similarity_k=15,
            sketch_size=100,
            seed=1729,
        )
    )
    return output


@pytest.fixture
def bundle(tmp_path: Path, bundle_template: Path) -> Path:
    output = tmp_path / "bundle"
    shutil.copytree(bundle_template, output)
    return output


@pytest.fixture
def circular_bundle(tmp_path: Path) -> Path:
    virus = tmp_path / "circular-viruses.fna"
    host = tmp_path / "circular-hosts.fna"
    metadata = tmp_path / "circular-metadata.tsv"
    virus.write_text(f">cv1\n{_dna(11, 31)}\n>cv2\n{_dna(12, 31)}\n", encoding="utf-8")
    host.write_text(f">ch1\n{_dna(13, 31)}\n>ch2\n{_dna(14, 31)}\n", encoding="utf-8")
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\ttopology\n"
        "cv1\tcv1\tvirus\tcircular\n"
        "cv2\tcv2\tvirus\tcircular\n"
        "ch1\tch1\thost\tcircular\n"
        "ch2\tch2\thost\tcircular\n",
        encoding="utf-8",
    )
    output = tmp_path / "circular-bundle"
    generate_benchmark(
        BenchmarkConfig(
            virus_paths=(virus,),
            host_paths=(host,),
            metadata_path=metadata,
            output_dir=output,
            splits=(SplitKind.RANDOM,),
            fragment_lengths=(31,),
            fragments_per_genome=2,
            seed=2718,
        )
    )
    return output


def _read_tsv(path: Path) -> tuple[list[dict[str, str]], tuple[str, ...]]:
    opener = gzip.open if path.name.endswith(".gz") else Path.open
    if path.name.endswith(".gz"):
        handle = opener(path, "rt", encoding="utf-8", newline="")
    else:
        handle = opener(path, "r", encoding="utf-8", newline="")
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), tuple(reader.fieldnames or ())


def _rewrite_tsv(
    path: Path,
    mutate: Callable[[list[dict[str, str]]], None],
) -> None:
    rows, columns = _read_tsv(path)
    mutate(rows)
    write_tsv(path, rows, columns)


def _rewrite_json(path: Path, mutate: Callable[[dict[str, Any]], None]) -> None:
    document = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(document, dict)
    mutate(document)
    write_json(path, document)


def _rewrite_fasta(path: Path, mutate: Callable[[list[Contig]], None]) -> None:
    records = list(read_fasta(path))
    mutate(records)
    write_fasta(records, path, overwrite=True)


def _sync_split_manifest(
    root: Path,
    kind: SplitKind,
    mutate: Callable[[dict[str, Any]], None],
) -> None:
    split_path = root / kind.directory_name / "split.json"
    split_manifest = json.loads(split_path.read_text(encoding="utf-8"))
    assert isinstance(split_manifest, dict)
    mutate(split_manifest)
    write_json(split_path, split_manifest)

    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["splits"][kind.value] = split_manifest
    write_json(manifest_path, manifest)


def _partition_statistics(
    records: list[Contig],
    rows: list[dict[str, str]],
) -> dict[str, object]:
    sequences = {record.sequence_id: record.sequence for record in records}
    lengths = Counter(len(sequence) for sequence in sequences.values())
    labels = Counter(row["label"] for row in rows)
    label_lengths = Counter(f"{row['label']}:{len(sequences[row['sequence_id']])}" for row in rows)
    genomes = Counter(row["source_genome_id"] for row in rows)
    total_bases = sum(map(len, sequences.values()))
    gc_bases = sum(sequence.count("G") + sequence.count("C") for sequence in sequences.values())
    ambiguous_bases = sum(
        sum(base not in "ACGT" for base in sequence) for sequence in sequences.values()
    )
    return {
        "records": len(rows),
        "bases": total_bases,
        "gc_fraction": gc_bases / total_bases,
        "ambiguous_fraction": ambiguous_bases / total_bases,
        "records_by_label": dict(sorted(labels.items())),
        "records_by_label_and_length": dict(sorted(label_lengths.items())),
        "records_by_length": {str(length): count for length, count in sorted(lengths.items())},
        "source_genomes": len(genomes),
        "records_by_genome": dict(sorted(genomes.items())),
    }


def _refresh_checksums(bundle: Path) -> None:
    write_checksums(bundle)


def test_valid_bundle_returns_structured_report(bundle: Path) -> None:
    report = validate_bundle(bundle)

    assert isinstance(report, ValidationReport)
    assert report.status == "pass"
    assert {kind for kind, _, _ in report.split_counts} == {kind.value for kind in SplitKind}
    primary_records = sum(train + test for _, train, test in report.split_counts)
    assert report.primary_fasta_records_verified == primary_records
    assert report.primary_truth_rows_verified == primary_records
    assert report.auxiliary_fasta_records_verified > 0
    assert report.auxiliary_truth_rows_verified == report.auxiliary_fasta_records_verified
    assert report.fasta_records_verified == report.truth_rows_verified
    assert report.fasta_records_verified == (
        report.primary_fasta_records_verified + report.auxiliary_fasta_records_verified
    )
    assert report.checksums_verified > 0
    serialized = report.as_dict()
    assert serialized["primary_fasta_records_verified"] == primary_records
    assert serialized["auxiliary_fasta_records_verified"] == (
        report.auxiliary_fasta_records_verified
    )
    assert serialized["splits"]["temporal"]["train_records"] > 0  # type: ignore[index]
    assert "Validated CHIMERA bundle" in report.summary()
    assert str(report) == report.summary()


def test_checksum_corruption_is_rejected_before_content_validation(bundle: Path) -> None:
    report = bundle / "REPORT.md"
    report.write_text(report.read_text(encoding="utf-8") + "tampered\n", encoding="utf-8")

    with pytest.raises(IntegrityError, match="Checksum verification failed"):
        validate_bundle(bundle)


def test_missing_embedded_schema_is_rejected_after_checksum_refresh(bundle: Path) -> None:
    schema = next((bundle / "schemas").glob("*.schema.json"))
    schema.unlink()
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="embedded schema inventory is incomplete"):
        validate_bundle(bundle)


def test_changed_embedded_schema_is_rejected_after_checksum_refresh(bundle: Path) -> None:
    schema = bundle / "schemas" / "truth-row.schema.json"

    def corrupt(document: dict[str, Any]) -> None:
        document["title"] = "attacker-controlled truth schema"

    _rewrite_json(schema, corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="embedded schema differs"):
        validate_bundle(bundle)


def test_root_manifest_missing_schema_required_field_is_rejected(bundle: Path) -> None:
    def corrupt(document: dict[str, Any]) -> None:
        document.pop("data_model")

    _rewrite_json(bundle / "manifest.json", corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="bundle schema violation"):
        validate_bundle(bundle)


def test_bad_bundle_marker_is_rejected_after_checksum_refresh(bundle: Path) -> None:
    (bundle / ".chimera-bundle").write_text("urn:chimera:benchmark-bundle:999\n", encoding="ascii")
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="marker must contain exactly"):
        validate_bundle(bundle)


def test_signed_extra_root_artifact_is_rejected_as_noncanonical(bundle: Path) -> None:
    (bundle / "attacker-controlled.txt").write_text("signed but undeclared\n", encoding="utf-8")
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="non-canonical bundle layout"):
        validate_bundle(bundle)


@pytest.mark.parametrize(
    ("field", "replacement", "message"),
    [
        ("release_date", "2099-01-01", "Genome-level release_date"),
        ("topology", "circular", "source sequence digest is incorrect"),
        ("taxonomy", '{"family":"Forged"}', "Genome taxonomy"),
    ],
)
def test_false_per_sequence_semantics_are_rejected(
    bundle: Path,
    field: str,
    replacement: str,
    message: str,
) -> None:
    def corrupt(rows: list[dict[str, str]]) -> None:
        rows[0][field] = replacement

    _rewrite_tsv(bundle / "sequences.tsv", corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match=message):
        validate_bundle(bundle)


def test_truth_valid_length_but_wrong_source_slice_is_rejected(bundle: Path) -> None:
    truth_path = bundle / SplitKind.RANDOM.directory_name / "train.truth.tsv.gz"
    sequence_rows, _ = _read_tsv(bundle / "sequences.tsv")
    lengths = {row["sequence_id"]: int(row["length_nt"]) for row in sequence_rows}

    def corrupt(rows: list[dict[str, str]]) -> None:
        row = next(
            item for item in rows if int(item["source_end"]) < lengths[item["source_sequence_id"]]
        )
        row["source_start"] = str(int(row["source_start"]) + 1)
        row["source_end"] = str(int(row["source_end"]) + 1)

    _rewrite_tsv(truth_path, corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="cannot be reconstructed"):
        validate_bundle(bundle)


def test_truth_reverse_complement_claim_is_reconstructed(bundle: Path) -> None:
    truth_path = bundle / SplitKind.RANDOM.directory_name / "train.truth.tsv.gz"

    def corrupt(rows: list[dict[str, str]]) -> None:
        rows[0]["strand"] = "+" if rows[0]["strand"] == "-" else "-"

    _rewrite_tsv(truth_path, corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="cannot be reconstructed"):
        validate_bundle(bundle)


def test_truth_circular_slice_is_reconstructed(circular_bundle: Path) -> None:
    truth_path = circular_bundle / SplitKind.RANDOM.directory_name / "train.truth.tsv.gz"

    def corrupt(rows: list[dict[str, str]]) -> None:
        row = rows[0]
        assert row["coordinate_system"] == "0-based-half-open-circular"
        source_length = 31
        start = (int(row["source_start"]) + 1) % source_length
        row["source_start"] = str(start)
        row["source_end"] = str(start + int(row["fragment_length"]))

    _rewrite_tsv(truth_path, corrupt)
    _refresh_checksums(circular_bundle)

    with pytest.raises(IntegrityError, match="cannot be reconstructed"):
        validate_bundle(circular_bundle)


@pytest.mark.parametrize("artifact", ["fasta", "truth"])
def test_fasta_truth_order_is_authenticated(bundle: Path, artifact: str) -> None:
    split_dir = bundle / SplitKind.GENOME.directory_name
    if artifact == "fasta":

        def reverse_records(records: list[Contig]) -> None:
            records.reverse()

        _rewrite_fasta(split_dir / "train.fasta.gz", reverse_records)
    else:

        def reverse_rows(rows: list[dict[str, str]]) -> None:
            rows.reverse()

        _rewrite_tsv(split_dir / "train.truth.tsv.gz", reverse_rows)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="identical ordered IDs"):
        validate_bundle(bundle)


def test_missing_assigned_genome_is_rejected_even_with_adjusted_statistics(bundle: Path) -> None:
    kind = SplitKind.RANDOM
    split_dir = bundle / kind.directory_name
    truth_path = split_dir / "train.truth.tsv.gz"
    fasta_path = split_dir / "train.fasta.gz"
    rows, columns = _read_tsv(truth_path)
    target_genome = rows[0]["source_genome_id"]
    retained_rows = [row for row in rows if row["source_genome_id"] != target_genome]
    retained_ids = {row["sequence_id"] for row in retained_rows}
    retained_records = [
        record for record in read_fasta(fasta_path) if record.sequence_id in retained_ids
    ]
    write_tsv(truth_path, retained_rows, columns)
    write_fasta(retained_records, fasta_path, overwrite=True)
    adjusted = _partition_statistics(retained_records, retained_rows)

    def adjust_manifest(document: dict[str, Any]) -> None:
        document["train"] = adjusted
        document["truth_rows"]["train"] = len(retained_rows)

    _sync_split_manifest(bundle, kind, adjust_manifest)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="emitted fragment IDs/order disagree"):
        validate_bundle(bundle)


@pytest.mark.parametrize(
    "statistic",
    ["gc_fraction", "ambiguous_fraction", "records_by_label_and_length"],
)
def test_composition_statistics_are_independently_recomputed(
    bundle: Path,
    statistic: str,
) -> None:
    def corrupt(document: dict[str, Any]) -> None:
        train = document["train"]
        if statistic == "records_by_label_and_length":
            key = next(iter(train[statistic]))
            train[statistic][key] += 1
        else:
            train[statistic] = (float(train[statistic]) + 0.25) % 1.0

    _sync_split_manifest(bundle, SplitKind.GENOME, corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match=rf"recorded train\.{statistic}"):
        validate_bundle(bundle)


def test_same_length_source_substitution_cannot_be_hidden_by_refreshing_digest_chain(
    circular_bundle: Path,
) -> None:
    truth_rows, _ = _read_tsv(
        circular_bundle / SplitKind.RANDOM.directory_name / "train.truth.tsv.gz"
    )
    target_truth = truth_rows[0]
    source_id = target_truth["source_sequence_id"]
    genome_id = target_truth["source_genome_id"]
    source_rows, _ = _read_tsv(circular_bundle / "sequences.tsv")
    target_source = next(row for row in source_rows if row["sequence_id"] == source_id)
    source_receipt = target_source["source_input_id"]
    source_records = list(read_fasta(circular_bundle / "source-sequences.fasta.gz"))
    original = next(record for record in source_records if record.sequence_id == source_id)
    position = int(target_truth["source_start"]) % len(original.sequence)
    replacement = "C" if original.sequence[position] != "C" else "A"
    mutated_sequence = (
        original.sequence[:position] + replacement + original.sequence[position + 1 :]
    )

    def replace_target(records: list[Contig]) -> None:
        index = next(i for i, record in enumerate(records) if record.sequence_id == source_id)
        records[index] = Contig(source_id, mutated_sequence)

    _rewrite_fasta(circular_bundle / "source-sequences.fasta.gz", replace_target)
    external_source = circular_bundle.parent / (
        "circular-viruses.fna" if target_source["label"] == "virus" else "circular-hosts.fna"
    )
    _rewrite_fasta(external_source, replace_target)

    exact_digest = hashlib.sha256(mutated_sequence.encode("ascii")).hexdigest()
    topology = target_source["topology"]
    canonical_digest = canonical_sequence_hash(
        mutated_sequence,
        circular=topology == "circular",
    )
    mutated_contig = Contig(source_id, mutated_sequence, topology=cast(Topology, topology))
    genome_digest = deterministic_genome_hash((mutated_contig,))
    external_digest = hashlib.sha256(external_source.read_bytes()).hexdigest()
    replacement_receipt = f"sha256:{external_digest}"

    def update_sequence_inventory(rows: list[dict[str, str]]) -> None:
        row = next(item for item in rows if item["sequence_id"] == source_id)
        row["sha256"] = exact_digest
        row["canonical_sha256"] = canonical_digest
        for item in rows:
            if item["source_input_id"] == source_receipt:
                item["source_input_id"] = replacement_receipt

    def update_reference_inventory(rows: list[dict[str, str]]) -> None:
        row = next(item for item in rows if item["genome_id"] == genome_id)
        assert row["contig_count"] == "1"
        row["sha256"] = genome_digest
        for item in rows:
            receipts = json.loads(item["source_input_ids"])
            item["source_input_ids"] = json.dumps(
                [replacement_receipt if value == source_receipt else value for value in receipts],
                separators=(",", ":"),
            )

    _rewrite_tsv(circular_bundle / "sequences.tsv", update_sequence_inventory)
    _rewrite_tsv(circular_bundle / "references.tsv", update_reference_inventory)

    split_dir = circular_bundle / SplitKind.RANDOM.directory_name

    def update_assignment(rows: list[dict[str, str]]) -> None:
        row = next(item for item in rows if item["genome_id"] == genome_id)
        row["group_id"] = f"sha256:{genome_digest}"

    _rewrite_tsv(split_dir / "assignments.tsv", update_assignment)
    for truth_path in split_dir.glob("*.truth.tsv.gz"):

        def update_truth(rows: list[dict[str, str]]) -> None:
            for row in rows:
                if row["source_genome_id"] == genome_id:
                    row["source_content_group_id"] = f"sha256:{genome_digest}"

        _rewrite_tsv(truth_path, update_truth)

    def update_manifest(document: dict[str, Any]) -> None:
        item = next(
            entry
            for entry in document["references"]["inputs"]
            if entry["content_id"] == source_receipt
        )
        item["content_id"] = replacement_receipt
        item["sha256"] = external_digest

    _rewrite_json(circular_bundle / "manifest.json", update_manifest)

    def update_resolved_config(document: dict[str, Any]) -> None:
        key = "virus_paths" if target_source["label"] == "virus" else "host_paths"
        document[key] = [
            replacement_receipt if value == source_receipt else value for value in document[key]
        ]

    _rewrite_json(circular_bundle / "resolved-config.json", update_resolved_config)
    _refresh_checksums(circular_bundle)

    with pytest.raises(IntegrityError, match="cannot be reconstructed"):
        validate_bundle(circular_bundle)


def test_root_exclusion_count_is_authenticated(bundle: Path) -> None:
    def corrupt(document: dict[str, Any]) -> None:
        document["references"]["preflight_exclusions"] += 1

    _rewrite_json(bundle / "manifest.json", corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="preflight_exclusions disagrees"):
        validate_bundle(bundle)


def test_root_exclusion_fields_are_authenticated(bundle: Path) -> None:
    references, _ = _read_tsv(bundle / "references.tsv")
    representative = references[0]
    rows, columns = _read_tsv(bundle / "excluded.tsv")
    forged = dict.fromkeys(columns, "")
    forged.update(
        {
            "genome_id": "forged-duplicate",
            "label": representative["label"],
            "split": "reference_preflight",
            "reason": "forged_exclusion_reason",
            "duplicate_of": representative["genome_id"],
            "source_sha256": representative["sha256"],
        }
    )
    rows.append(forged)
    write_tsv(bundle / "excluded.tsv", rows, columns)

    def adjust_count(document: dict[str, Any]) -> None:
        document["references"]["preflight_exclusions"] = len(rows)

    _rewrite_json(bundle / "manifest.json", adjust_count)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="invalid preflight duplicate exclusion provenance"):
        validate_bundle(bundle)


@pytest.mark.parametrize(
    ("drift", "message"),
    [
        ("candidate_label", "label disagrees with source inventory"),
        ("candidate_source", "source identifiers do not resolve"),
        ("stratum_row", "similarity truth field"),
    ],
)
def test_similarity_auxiliary_views_reject_label_source_and_row_drift(
    bundle: Path,
    drift: str,
    message: str,
) -> None:
    split_dir = bundle / SplitKind.SIMILARITY.directory_name
    candidate_path = split_dir / "candidate_test.truth.tsv.gz"
    candidate_rows, _ = _read_tsv(candidate_path)
    if drift == "candidate_label":

        def corrupt_candidate_label(rows: list[dict[str, str]]) -> None:
            rows[0]["label"] = "host" if rows[0]["label"] == "virus" else "virus"

        _rewrite_tsv(candidate_path, corrupt_candidate_label)
    elif drift == "candidate_source":
        references, _ = _read_tsv(bundle / "references.tsv")
        target = candidate_rows[0]
        replacement = next(
            row["genome_id"]
            for row in references
            if row["label"] == target["label"] and row["genome_id"] != target["source_genome_id"]
        )

        def corrupt_candidate_source(rows: list[dict[str, str]]) -> None:
            rows[0]["source_genome_id"] = replacement

        _rewrite_tsv(candidate_path, corrupt_candidate_source)
    else:
        occupied_bin = candidate_rows[0]["similarity_bin"]
        stratum_path = split_dir / "test_strata" / f"{occupied_bin}.truth.tsv.gz"

        def corrupt_stratum_row(rows: list[dict[str, str]]) -> None:
            rows[0]["nearest_train_genome_id"] = "forged-neighbor"

        _rewrite_tsv(stratum_path, corrupt_stratum_row)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match=message):
        validate_bundle(bundle)


@pytest.mark.parametrize(
    ("drift", "message"), [("seed", "seed disagrees"), ("splits", "split set disagrees")]
)
def test_resolved_config_seed_and_split_set_are_authenticated(
    bundle: Path,
    drift: str,
    message: str,
) -> None:
    def corrupt(document: dict[str, Any]) -> None:
        if drift == "seed":
            document["seed"] += 1
        else:
            document["splits"].remove(SplitKind.TAXONOMY.value)

    _rewrite_json(bundle / "resolved-config.json", corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match=message):
        validate_bundle(bundle)


@pytest.mark.parametrize(
    ("corruption", "message"),
    [
        ("unexpected_field", "incomplete or has unexpected fields"),
        ("naive_timestamp", "must include a UTC offset"),
        ("non_string_command", "command must be a JSON array of strings"),
    ],
)
def test_malformed_execution_record_is_rejected(
    bundle: Path,
    corruption: str,
    message: str,
) -> None:
    def corrupt(document: dict[str, Any]) -> None:
        if corruption == "unexpected_field":
            document["attacker_note"] = "trusted"
        elif corruption == "naive_timestamp":
            document["started_at_utc"] = "2026-08-08T12:00:00"
        else:
            document["command"] = ["chimera", 7]

    _rewrite_json(bundle / "execution.json", corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match=message):
        validate_bundle(bundle)


def test_fasta_truth_identifier_drift_is_rejected_even_with_fresh_checksums(
    bundle: Path,
) -> None:
    truth = bundle / SplitKind.GENOME.directory_name / "test.truth.tsv.gz"

    def corrupt(rows: list[dict[str, str]]) -> None:
        rows[0]["sequence_id"] = "frag-00000000000000000000000000000000"

    _rewrite_tsv(truth, corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="FASTA and truth must have identical ordered IDs"):
        validate_bundle(bundle)


def test_truth_coordinate_and_sequence_length_corruption_is_rejected(bundle: Path) -> None:
    truth = bundle / SplitKind.RANDOM.directory_name / "train.truth.tsv.gz"

    def corrupt(rows: list[dict[str, str]]) -> None:
        rows[0]["source_end"] = str(int(rows[0]["source_end"]) + 1)

    _rewrite_tsv(truth, corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="half-open coordinates"):
        validate_bundle(bundle)


def test_truth_source_genome_must_resolve_with_source_sequence(bundle: Path) -> None:
    split_dir = bundle / SplitKind.GENOME.directory_name

    def corrupt(rows: list[dict[str, str]]) -> None:
        target_genome = rows[0]["source_genome_id"]
        for row in rows:
            if row["source_genome_id"] == target_genome:
                row["source_genome_id"] = "missing-source"

    _rewrite_tsv(split_dir / "test.truth.tsv.gz", corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="source identifiers do not resolve"):
        validate_bundle(bundle)


def test_assignment_release_date_must_match_reference_inventory(bundle: Path) -> None:
    assignments = bundle / SplitKind.TEMPORAL.directory_name / "assignments.tsv"

    def corrupt(rows: list[dict[str, str]]) -> None:
        training = next(row for row in rows if row["partition"] == "train")
        training["release_date"] = "2099-01-01"

    _rewrite_tsv(assignments, corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match=r"release_date disagrees with references.tsv"):
        validate_bundle(bundle)


def test_reference_taxonomy_rejects_legacy_flat_serialization(bundle: Path) -> None:
    references_path = bundle / "references.tsv"

    def corrupt_references(rows: list[dict[str, str]]) -> None:
        rows[0]["taxonomy"] = "family=Alpha"

    _rewrite_tsv(references_path, corrupt_references)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="taxonomy must be a JSON object"):
        validate_bundle(bundle)


def test_similarity_evidence_is_recomputed_from_resolved_inputs(bundle: Path) -> None:
    split_dir = bundle / SplitKind.SIMILARITY.directory_name

    def corrupt(rows: list[dict[str, str]]) -> None:
        retained = next(row for row in rows if row["partition"] == "test")
        retained["similarity_method"] = "tampered-method"

    _rewrite_tsv(split_dir / "assignments.tsv", corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(
        IntegrityError,
        match="similarity_method disagrees with independently recomputed split assignment",
    ):
        validate_bundle(bundle)


def test_similarity_candidate_truth_must_match_assignments(bundle: Path) -> None:
    split_dir = bundle / SplitKind.SIMILARITY.directory_name

    def corrupt(rows: list[dict[str, str]]) -> None:
        rows[0]["similarity_bin"] = (
            "low_similarity"
            if rows[0]["similarity_bin"] != "low_similarity"
            else "moderate_similarity"
        )

    _rewrite_tsv(split_dir / "candidate_test.truth.tsv.gz", corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="similarity truth field"):
        validate_bundle(bundle)


def test_similarity_strata_must_partition_candidate_test(bundle: Path) -> None:
    split_dir = bundle / SplitKind.SIMILARITY.directory_name
    candidate_rows, _ = _read_tsv(split_dir / "candidate_test.truth.tsv.gz")
    occupied_bin = candidate_rows[0]["similarity_bin"]
    truth_path = split_dir / "test_strata" / f"{occupied_bin}.truth.tsv.gz"

    def corrupt(rows: list[dict[str, str]]) -> None:
        rows.pop()

    _rewrite_tsv(truth_path, corrupt)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="FASTA and truth must have identical ordered IDs"):
        validate_bundle(bundle)


def test_exclusion_row_must_be_backed_by_an_excluded_assignment(bundle: Path) -> None:
    split_dir = bundle / SplitKind.SIMILARITY.directory_name
    assignments, _ = _read_tsv(split_dir / "assignments.tsv")
    target = next(row for row in assignments if row["candidate_partition"] == "test")
    references, _ = _read_tsv(bundle / "references.tsv")
    reference = next(row for row in references if row["genome_id"] == target["genome_id"])
    excluded_path = split_dir / "excluded.tsv"
    _, columns = _read_tsv(excluded_path)
    fabricated = dict.fromkeys(columns, "")
    fabricated.update(
        {
            "genome_id": target["genome_id"],
            "label": target["label"],
            "split": "similarity",
            "reason": "fabricated_exclusion",
            "source_sha256": reference["sha256"],
            "source_accession_version": reference["accession_version"],
            "release_date": reference["release_date"],
        }
    )
    write_tsv(excluded_path, (fabricated,), columns)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match="not backed by an excluded assignment"):
        validate_bundle(bundle)


def test_recomputed_counts_override_recorded_pass_flag(bundle: Path) -> None:
    split_path = bundle / SplitKind.GENOME.directory_name / "split.json"
    root_path = bundle / "manifest.json"
    split_manifest = json.loads(split_path.read_text(encoding="utf-8"))
    root_manifest = json.loads(root_path.read_text(encoding="utf-8"))
    split_manifest["validation"]["status"] = "pass"
    split_manifest["train"]["records"] += 1
    root_manifest["splits"]["genome"] = split_manifest
    write_json(split_path, split_manifest)
    write_json(root_path, root_manifest)
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match=r"recorded train\.records"):
        validate_bundle(bundle)


def test_recorded_failure_flag_is_rejected_by_manifest_schema(bundle: Path) -> None:
    split_path = bundle / SplitKind.GENOME.directory_name / "split.json"
    root_path = bundle / "manifest.json"
    split_manifest = json.loads(split_path.read_text(encoding="utf-8"))
    root_manifest = json.loads(root_path.read_text(encoding="utf-8"))
    split_manifest["validation"]["status"] = "fail"
    root_manifest["splits"]["genome"] = split_manifest
    write_json(split_path, split_manifest)
    write_json(root_path, root_manifest)
    _refresh_checksums(bundle)

    with pytest.raises(SchemaValidationError, match="'pass' was expected"):
        validate_bundle(bundle)


def test_missing_declared_split_file_is_rejected_after_checksum_refresh(bundle: Path) -> None:
    missing = bundle / SplitKind.GENOME.directory_name / "test.truth.tsv.gz"
    missing.unlink()
    _refresh_checksums(bundle)

    with pytest.raises(IntegrityError, match=r"non-canonical bundle layout.*test\.truth\.tsv\.gz"):
        validate_bundle(bundle)


def test_validation_rejects_symlink_and_invalid_roots(bundle: Path, tmp_path: Path) -> None:
    link = tmp_path / "bundle-link"
    link.symlink_to(bundle, target_is_directory=True)
    with pytest.raises(IntegrityError, match="must not be a symbolic link"):
        validate_bundle(link)

    with pytest.raises(IntegrityError, match="Cannot access bundle root"):
        validate_bundle(tmp_path / "absent")

    regular_file = tmp_path / "not-a-directory"
    regular_file.write_text("not a bundle\n", encoding="utf-8")
    with pytest.raises(IntegrityError, match="not a directory"):
        validate_bundle(regular_file)


@pytest.mark.parametrize(
    ("value", "message"),
    [("not-an-integer", "must be an integer"), ("01", "canonical integer syntax")],
)
def test_validation_integer_parser_rejects_noncanonical_values(
    tmp_path: Path,
    value: str,
    message: str,
) -> None:
    with pytest.raises(IntegrityError, match=message):
        validation_module._parse_int(
            value,
            path=tmp_path / "row.tsv",
            line_number=2,
            field="count",
        )


@pytest.mark.parametrize(
    ("value", "allow_empty", "message"),
    [
        ("not-a-number", False, "must be a number"),
        ("nan", False, "must be finite"),
        ("1.01", False, r"on \[0, 1\]"),
        ("", False, "must be a number"),
    ],
)
def test_validation_fraction_parser_rejects_invalid_values(
    tmp_path: Path,
    value: str,
    allow_empty: bool,
    message: str,
) -> None:
    with pytest.raises(IntegrityError, match=message):
        validation_module._parse_fraction(
            value,
            path=tmp_path / "row.tsv",
            line_number=2,
            field="similarity",
            allow_empty=allow_empty,
        )


def test_validation_date_parser_rejects_non_iso_date(tmp_path: Path) -> None:
    with pytest.raises(IntegrityError, match="ISO YYYY-MM-DD"):
        validation_module._parse_date(
            "08/08/2026",
            path=tmp_path / "row.tsv",
            line_number=2,
            field="release_date",
        )


@pytest.mark.parametrize(
    ("value", "message"),
    [
        ("{", "must be a JSON object"),
        ('{"rank":"one","rank":"two"}', "duplicate taxonomy key"),
        ("[]", "must map strings to strings"),
        ('{"rank":7}', "must map strings to strings"),
    ],
)
def test_validation_json_map_parser_rejects_ambiguous_values(
    tmp_path: Path,
    value: str,
    message: str,
) -> None:
    with pytest.raises(IntegrityError, match=message):
        validation_module._parse_json_string_map(
            value,
            path=tmp_path / "row.tsv",
            line_number=2,
            field="taxonomy",
        )


@pytest.mark.parametrize(
    ("value", "message"),
    [
        ("[", "must be a JSON array"),
        ("[]", "nonempty array of unique strings"),
        ('["one","one"]', "nonempty array of unique strings"),
        ("[1]", "nonempty array of unique strings"),
    ],
)
def test_validation_json_list_parser_rejects_ambiguous_values(
    tmp_path: Path,
    value: str,
    message: str,
) -> None:
    with pytest.raises(IntegrityError, match=message):
        validation_module._parse_json_string_list(
            value,
            path=tmp_path / "row.tsv",
            line_number=2,
            field="sequence_ids",
        )


def test_validation_json_reader_rejects_malformed_documents(tmp_path: Path) -> None:
    cases = (
        ('{"key":1,"key":2}', "duplicate JSON key"),
        ("{", "Cannot read JSON object"),
        ("[]", "top-level JSON value must be an object"),
    )
    for index, (content, message) in enumerate(cases):
        path = tmp_path / f"case-{index}.json"
        path.write_text(content, encoding="utf-8")
        with pytest.raises(IntegrityError, match=message):
            validation_module._read_json(path)

    with pytest.raises(IntegrityError, match="Cannot read JSON object"):
        validation_module._read_json(tmp_path / "missing.json")


def test_validation_tsv_reader_rejects_structural_corruption(tmp_path: Path) -> None:
    cases = (
        ("required\trequired\nvalue\tvalue\n", "duplicate columns"),
        ("other\nvalue\n", "missing required column"),
        ("required\nvalue\textra\n", "more fields than the header"),
    )
    for index, (content, message) in enumerate(cases):
        path = tmp_path / f"case-{index}.tsv"
        path.write_text(content, encoding="utf-8")
        with pytest.raises(IntegrityError, match=message):
            validation_module._read_tsv(path, ("required",))

    undecodable = tmp_path / "invalid-utf8.tsv"
    undecodable.write_bytes(b"required\n\xff\n")
    with pytest.raises(IntegrityError, match="Cannot read TSV"):
        validation_module._read_tsv(undecodable, ("required",))


def test_validation_execution_rejects_remaining_invalid_fields(tmp_path: Path) -> None:
    manifest = {"randomness": {"python_version": "3.12.0"}}
    valid = {
        "started_at_utc": "2026-08-08T12:00:00+00:00",
        "finished_at_utc": "2026-08-08T12:01:00+00:00",
        "python": "3.12.0",
        "platform": "test-platform",
        "command": ["chimera"],
        "status": "complete",
    }
    cases = (
        ("timestamp_type", "must be an ISO timestamp"),
        ("timestamp_syntax", "is not an ISO timestamp"),
        ("reverse_time", "precedes started_at_utc"),
        ("python", "Python version disagrees"),
        ("platform", "platform must be a non-empty string"),
    )
    for case, message in cases:
        execution = json.loads(json.dumps(valid))
        if case == "timestamp_type":
            execution["started_at_utc"] = 7
        elif case == "timestamp_syntax":
            execution["started_at_utc"] = "not-a-timestamp"
        elif case == "reverse_time":
            execution["finished_at_utc"] = "2026-08-08T11:59:00+00:00"
        elif case == "python":
            execution["python"] = "0.0.0"
        else:
            execution["platform"] = ""
        write_json(tmp_path / "execution.json", execution)
        with pytest.raises(IntegrityError, match=message):
            validation_module._validate_execution(tmp_path, manifest)


def _valid_manifest_provenance() -> dict[str, Any]:
    digest_a = "a" * 64
    digest_b = "b" * 64
    return {
        "tool": {
            "name": "CHIMERA",
            "version": validation_module.__version__,
            "software_content_sha256": validation_module.software_content_sha256(),
            "git_revision": "unknown",
            "git_dirty": None,
        },
        "randomness": {
            "algorithm": "Python random.Random with semantic BLAKE2b-derived sub-seeds",
            "seed_derivation": "chimera.seed.v1",
            "python_implementation": "CPython",
            "python_version": "3.12.0",
        },
        "references": {
            "inputs": [
                {
                    "content_id": f"sha256:{digest_a}",
                    "role": "reference_fasta",
                    "sha256": digest_a,
                },
                {
                    "content_id": f"sha256:{digest_b}",
                    "role": "reference_fasta",
                    "sha256": digest_b,
                },
            ]
        },
    }


def test_validation_manifest_provenance_rejects_every_untrusted_component(
    tmp_path: Path,
) -> None:
    cases = (
        ("tool", "invalid tool provenance"),
        ("version", "tool version is unsupported"),
        ("software", "software content receipt disagrees"),
        ("revision", "invalid git_revision"),
        ("dirty_type", "git_dirty must be boolean or null"),
        ("unknown_dirty", "git_dirty must be null"),
        ("randomness", "invalid randomness provenance"),
        ("missing_inputs", "must inventory source inputs"),
        ("malformed_input", r"inputs\[0\] is malformed"),
        ("invalid_input", r"inputs\[0\] has invalid provenance"),
        ("duplicate_input", r"inputs\[1\] has invalid provenance"),
    )
    for case, message in cases:
        manifest = json.loads(json.dumps(_valid_manifest_provenance()))
        if case == "tool":
            manifest["tool"] = []
        elif case == "version":
            manifest["tool"]["version"] = "0.0.0"
        elif case == "software":
            manifest["tool"]["software_content_sha256"] = "0" * 64
        elif case == "revision":
            manifest["tool"]["git_revision"] = "not-a-revision"
        elif case == "dirty_type":
            manifest["tool"]["git_revision"] = "1" * 40
            manifest["tool"]["git_dirty"] = "yes"
        elif case == "unknown_dirty":
            manifest["tool"]["git_dirty"] = False
        elif case == "randomness":
            manifest["randomness"]["seed_derivation"] = "untrusted"
        elif case == "missing_inputs":
            manifest["references"]["inputs"] = []
        elif case == "malformed_input":
            manifest["references"]["inputs"][0]["extra"] = "field"
        elif case == "invalid_input":
            manifest["references"]["inputs"][0]["content_id"] = "sha256:wrong"
        else:
            manifest["references"]["inputs"][1] = dict(manifest["references"]["inputs"][0])
        with pytest.raises(IntegrityError, match=message):
            validation_module._validate_manifest_provenance(tmp_path, manifest)


def _source_for_input(input_id: str, label: str) -> Any:
    digest = "0" * 64
    return validation_module._SourceSequence(
        sequence_id=f"{label}-sequence",
        genome_id=f"{label}-genome",
        label=label,
        accession_version="",
        release_date="",
        topology="linear",
        sequence="A",
        sha256=digest,
        canonical_sha256=digest,
        source_input_id=input_id,
        taxonomy={},
        metadata_extra={},
    )


def _input_receipt(content_id: str, role: str = "reference_fasta") -> dict[str, str]:
    return {
        "content_id": content_id,
        "role": role,
        "sha256": content_id.removeprefix("sha256:"),
    }


def test_validation_input_receipt_links_reject_every_mismatch(tmp_path: Path) -> None:
    virus_id = f"sha256:{'a' * 64}"
    host_id = f"sha256:{'b' * 64}"
    extra_id = f"sha256:{'c' * 64}"
    metadata_id = f"sha256:{'d' * 64}"
    similarity_id = f"sha256:{'e' * 64}"

    def config(**changes: Any) -> BenchmarkConfig:
        values: dict[str, Any] = {
            "virus_paths": (Path(virus_id),),
            "host_paths": (Path(host_id),),
            "output_dir": tmp_path / "out",
        }
        values.update(changes)
        return BenchmarkConfig(**values)

    base_manifest = {"references": {"inputs": [_input_receipt(virus_id), _input_receipt(host_id)]}}
    base_sources = {
        "virus": _source_for_input(virus_id, "virus"),
        "host": _source_for_input(host_id, "host"),
    }
    cases = (
        ("virus", "viral input IDs disagree"),
        ("host", "host input IDs disagree"),
        ("overlap", "assigned to both biological labels"),
        ("receipts", "reference input receipts disagree"),
        ("metadata", "metadata input provenance disagrees"),
        ("similarity", "similarity input provenance disagrees"),
    )
    for case, message in cases:
        manifest = json.loads(json.dumps(base_manifest))
        sources = dict(base_sources)
        current_config = config()
        if case == "virus":
            sources["virus"] = _source_for_input(extra_id, "virus")
        elif case == "host":
            sources["host"] = _source_for_input(extra_id, "host")
        elif case == "overlap":
            current_config = config(host_paths=(Path(virus_id),))
            sources["host"] = _source_for_input(virus_id, "host")
            manifest["references"]["inputs"] = [_input_receipt(virus_id)]
        elif case == "receipts":
            manifest["references"]["inputs"].append(_input_receipt(extra_id))
        elif case == "metadata":
            current_config = config(metadata_path=Path(metadata_id))
        else:
            current_config = config(similarity_table=Path(similarity_id))
        with pytest.raises(IntegrityError, match=message):
            validation_module._validate_input_inventory_links(
                tmp_path,
                manifest,
                current_config,
                sources,
            )


def test_validation_exact_directory_wraps_enumeration_errors(tmp_path: Path) -> None:
    with pytest.raises(IntegrityError, match="Cannot enumerate bundle directory"):
        validation_module._validate_exact_directory(tmp_path / "missing", set())


def _reference(
    *,
    genome_id: str = "genome-a",
    label: str = "virus",
    digest: str | None = None,
    sequence_ids: frozenset[str] = frozenset({"sequence-a"}),
    source_input_ids: frozenset[str] = frozenset({f"sha256:{'b' * 64}"}),
    accession_version: str = "",
    release_date: str = "",
    taxonomy: dict[str, str] | None = None,
    contig_count: int = 1,
    length_nt: int = 1,
) -> Any:
    if digest is None:
        digest = deterministic_genome_hash((Contig(sequence_id="sequence-a", sequence="A"),))
    return validation_module._Reference(
        genome_id=genome_id,
        label=label,
        accession_version=accession_version,
        release_date=release_date,
        digest=digest,
        sequence_ids=sequence_ids,
        source_input_ids=source_input_ids,
        taxonomy=taxonomy or {},
        contig_count=contig_count,
        length_nt=length_nt,
        metadata_extra={},
    )


def _reference_row(**changes: str) -> dict[str, str]:
    row = {
        "genome_id": "genome-a",
        "label": "virus",
        "accession_version": "",
        "release_date": "",
        "sequence_ids": '["sequence-a"]',
        "contig_count": "1",
        "length_nt": "1",
        "sha256": deterministic_genome_hash((Contig(sequence_id="sequence-a", sequence="A"),)),
        "source_input_ids": f'["sha256:{"b" * 64}"]',
        "taxonomy": "{}",
        "metadata_extra": "{}",
    }
    row.update(changes)
    return row


def test_validation_reference_reader_defends_against_semantic_corruption(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(validation_module, "validate_instance", lambda *_args, **_kwargs: None)
    cases: tuple[tuple[list[dict[str, str]], object, str], ...] = (
        ([_reference_row(genome_id="")], {"count": 1}, "genome_id must not be empty"),
        ([_reference_row(), _reference_row()], {"count": 2}, "duplicate genome_id"),
        ([_reference_row(label="other")], {"count": 1}, "label must be 'virus' or 'host'"),
        ([_reference_row(sha256="bad")], {"count": 1}, "invalid canonical genome SHA-256"),
        (
            [_reference_row(), _reference_row(genome_id="genome-b")],
            {"count": 2},
            "duplicate canonical genome content",
        ),
        (
            [_reference_row(source_input_ids='["input.fna"]')],
            {"count": 1},
            "source_input_ids must be content-addressed",
        ),
        ([_reference_row(contig_count="0")], {"count": 1}, "must be positive"),
        (
            [_reference_row(contig_count="2")],
            {"count": 1},
            "contig_count disagrees with sequence_ids",
        ),
        ([_reference_row()], [], "references must be an object"),
        ([_reference_row()], {"count": 2}, "reference count does not match"),
    )
    path = tmp_path / "references.tsv"
    for rows, reference_manifest, message in cases:
        write_tsv(path, rows, validation_module.REFERENCE_COLUMNS)
        with pytest.raises(IntegrityError, match=message):
            validation_module._read_references(
                tmp_path,
                {"references": reference_manifest},
            )


def _assignment_row(reference: Any, **changes: str) -> dict[str, str]:
    row = dict.fromkeys(validation_module._ASSIGNMENT_COLUMNS, "")
    row.update(
        {
            "genome_id": reference.genome_id,
            "group_id": f"sha256:{reference.digest}",
            "label": reference.label,
            "partition": "train",
            "reason": "test",
            "release_date": reference.release_date,
        }
    )
    row.update(changes)
    return row


def test_validation_assignment_reader_defends_against_semantic_corruption(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(validation_module, "validate_instance", lambda *_args, **_kwargs: None)
    reference = _reference()
    references = {reference.genome_id: reference}
    base = _assignment_row(reference)
    cases: tuple[tuple[list[dict[str, str]], str], ...] = (
        ([base, dict(base)], "duplicate genome assignment"),
        ([_assignment_row(reference, genome_id="unknown")], "assignment names unknown genome"),
        ([_assignment_row(reference, partition="both")], "partition 'both' is invalid"),
        ([_assignment_row(reference, label="host")], "label disagrees with references.tsv"),
        (
            [_assignment_row(reference, group_id=f"sha256:{'0' * 64}")],
            "content group disagrees",
        ),
        ([_assignment_row(reference, release_date="2026-01-01")], "release_date disagrees"),
        ([_assignment_row(reference, similarity_bin="mystery")], "unknown similarity_bin"),
        (
            [_assignment_row(reference, max_train_similarity="0.10")],
            "canonical round-trip syntax",
        ),
        ([], "assignments must cover every reference exactly once"),
    )
    path = tmp_path / "assignments.tsv"
    for rows, message in cases:
        write_tsv(path, rows, validation_module._ASSIGNMENT_COLUMNS)
        with pytest.raises(IntegrityError, match=message):
            validation_module._read_assignments(path, SplitKind.GENOME, references)


def _exclusion_row(reference: Any, **changes: str) -> dict[str, str]:
    row = dict.fromkeys(validation_module._EXCLUSION_COLUMNS, "")
    row.update(
        {
            "genome_id": reference.genome_id,
            "label": reference.label,
            "split": SplitKind.GENOME.value,
            "reason": "excluded-for-test",
            "source_sha256": reference.digest,
            "source_accession_version": reference.accession_version,
            "release_date": reference.release_date,
        }
    )
    row.update(changes)
    return row


def test_validation_split_exclusion_reader_checks_all_linked_fields(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(validation_module, "validate_instance", lambda *_args, **_kwargs: None)
    reference = _reference()
    references = {reference.genome_id: reference}
    assignment = _assignment_row(
        reference,
        partition="excluded",
        reason="excluded-for-test",
    )
    assignments = {reference.genome_id: assignment}
    base = _exclusion_row(reference)
    cases: tuple[tuple[list[dict[str, str]], str], ...] = (
        ([base, dict(base)], "duplicate excluded genome"),
        (
            [_exclusion_row(reference, genome_id="unknown")],
            "not backed by an excluded assignment",
        ),
        ([_exclusion_row(reference, split="temporal")], "split/label disagrees"),
        ([_exclusion_row(reference, reason="wrong")], "reason disagrees"),
        (
            [_exclusion_row(reference, source_sha256="0" * 64)],
            "source_sha256 disagrees",
        ),
        ([], "account for every excluded assignment exactly once"),
    )
    path = tmp_path / "excluded.tsv"
    for rows, message in cases:
        write_tsv(path, rows, validation_module._EXCLUSION_COLUMNS)
        with pytest.raises(IntegrityError, match=message):
            validation_module._validate_exclusions(
                path,
                SplitKind.GENOME,
                assignments,
                references,
            )


def test_validation_partition_stats_require_objects_and_truth_counts(tmp_path: Path) -> None:
    partition = validation_module._Partition(
        name="train",
        rows=(),
        sequences={},
        ordered_ids=(),
    )
    with pytest.raises(IntegrityError, match="missing 'train' statistics object"):
        validation_module._compare_partition_stats(tmp_path / "split.json", {}, partition)

    empty_stats = {
        "records": 0,
        "bases": 0,
        "gc_fraction": None,
        "ambiguous_fraction": None,
        "records_by_label": {},
        "records_by_label_and_length": {},
        "records_by_length": {},
        "records_by_genome": {},
        "source_genomes": 0,
    }
    with pytest.raises(IntegrityError, match=r"recorded truth_rows\.train is incorrect"):
        validation_module._compare_partition_stats(
            tmp_path / "split.json",
            {"train": empty_stats},
            partition,
        )


def test_validation_execution_rejects_shape_timezone_and_command(tmp_path: Path) -> None:
    manifest = {"randomness": {"python_version": "3.12.0"}}
    base = {
        "started_at_utc": "2026-08-08T12:00:00+00:00",
        "finished_at_utc": "2026-08-08T12:01:00+00:00",
        "python": "3.12.0",
        "platform": "test-platform",
        "command": ["chimera"],
        "status": "complete",
    }
    cases = (
        ({**base, "status": "running"}, "execution record is incomplete"),
        ({**base, "started_at_utc": "2026-08-08T12:00:00"}, "must include a UTC offset"),
        ({**base, "command": [1]}, "command must be a JSON array of strings"),
    )
    for execution, message in cases:
        write_json(tmp_path / "execution.json", execution)
        with pytest.raises(IntegrityError, match=message):
            validation_module._validate_execution(tmp_path, manifest)


def _source_row(**changes: str) -> dict[str, str]:
    sequence = changes.pop("_sequence", "A")
    topology = changes.get("topology", "linear")
    row = {
        "sequence_id": "sequence-a",
        "genome_id": "genome-a",
        "label": "virus",
        "accession_version": "",
        "release_date": "",
        "topology": topology,
        "length_nt": str(len(sequence)),
        "sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        "canonical_sha256": canonical_sequence_hash(
            sequence,
            circular=topology == "circular",
        ),
        "source_input_id": f"sha256:{'b' * 64}",
        "taxonomy": "{}",
        "metadata_extra": "{}",
    }
    row.update(changes)
    return row


def test_validation_source_inventory_reader_checks_independent_semantics(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(validation_module, "validate_instance", lambda *_args, **_kwargs: None)
    table = tmp_path / "sequences.tsv"
    fasta = tmp_path / "source-sequences.fasta.gz"
    base_contig = Contig(sequence_id="sequence-a", sequence="A")
    base_reference = _reference()

    def validate(
        rows: list[dict[str, str]],
        contigs: list[Contig],
        references: dict[str, Any],
        message: str,
    ) -> None:
        write_tsv(table, rows, validation_module.SEQUENCE_COLUMNS)
        write_fasta(contigs, fasta, overwrite=True)
        with pytest.raises(IntegrityError, match=message):
            validation_module._read_source_sequences(tmp_path, references)

    with gzip.open(fasta, "wt", encoding="utf-8") as handle:
        handle.write("not FASTA\n")
    write_tsv(table, [_source_row()], validation_module.SEQUENCE_COLUMNS)
    with pytest.raises(IntegrityError, match="Cannot validate source FASTA"):
        validation_module._read_source_sequences(
            tmp_path,
            {base_reference.genome_id: base_reference},
        )

    cases = (
        (
            [_source_row()],
            [Contig(sequence_id="different", sequence="A")],
            {base_reference.genome_id: base_reference},
            "identical ordered IDs",
        ),
        ([_source_row()], [base_contig], {}, "unknown genome_id"),
        (
            [_source_row(label="host")],
            [base_contig],
            {base_reference.genome_id: base_reference},
            "label disagrees",
        ),
        (
            [_source_row()],
            [base_contig],
            {
                base_reference.genome_id: _reference(
                    sequence_ids=frozenset({"different"}),
                )
            },
            "sequence is absent from its reference group",
        ),
        (
            [_source_row(topology="unknown")],
            [base_contig],
            {base_reference.genome_id: base_reference},
            "topology must be linear or circular",
        ),
        (
            [_source_row(length_nt="2")],
            [base_contig],
            {base_reference.genome_id: base_reference},
            "length_nt disagrees",
        ),
        (
            [_source_row(source_input_id="input.fna")],
            [base_contig],
            {base_reference.genome_id: base_reference},
            "source_input_id must be content-addressed",
        ),
        (
            [_source_row()],
            [base_contig],
            {
                base_reference.genome_id: _reference(
                    sequence_ids=frozenset({"sequence-a", "missing"}),
                )
            },
            "inventory every retained source sequence",
        ),
        (
            [_source_row()],
            [base_contig],
            {base_reference.genome_id: _reference(contig_count=2)},
            "Source contig count disagrees",
        ),
        (
            [_source_row()],
            [base_contig],
            {base_reference.genome_id: _reference(length_nt=2)},
            "Source length disagrees",
        ),
        (
            [_source_row()],
            [base_contig],
            {base_reference.genome_id: _reference(digest="0" * 64)},
            "Source content digest disagrees",
        ),
        (
            [_source_row()],
            [base_contig],
            {
                base_reference.genome_id: _reference(
                    source_input_ids=frozenset({f"sha256:{'c' * 64}"})
                )
            },
            "Reference source input IDs disagree",
        ),
        (
            [_source_row()],
            [base_contig],
            {
                base_reference.genome_id: _reference(
                    accession_version="ACCESSION.1",
                )
            },
            "accession_version.*disagrees",
        ),
    )
    for rows, contigs, references, message in cases:
        validate(rows, contigs, references, message)

    conflicting_contigs = [
        Contig(
            sequence_id="sequence-a",
            sequence="A",
            taxonomy=(("family", "Alpha"),),
        ),
        Contig(
            sequence_id="sequence-b",
            sequence="C",
            taxonomy=(("family", "Beta"),),
        ),
    ]
    conflicting_rows = [
        _source_row(taxonomy='{"family":"Alpha"}'),
        _source_row(
            sequence_id="sequence-b",
            taxonomy='{"family":"Beta"}',
            _sequence="C",
        ),
    ]
    conflicting_reference = _reference(
        digest=deterministic_genome_hash(tuple(conflicting_contigs)),
        sequence_ids=frozenset({"sequence-a", "sequence-b"}),
        contig_count=2,
        length_nt=2,
    )
    validate(
        conflicting_rows,
        conflicting_contigs,
        {conflicting_reference.genome_id: conflicting_reference},
        "Per-sequence taxonomy conflicts",
    )


def _fragment_source(*, topology: str = "linear") -> Any:
    sequence = "ACGT"
    return validation_module._SourceSequence(
        sequence_id="sequence-a",
        genome_id="genome-a",
        label="virus",
        accession_version="ACCESSION.1",
        release_date="2026-01-01",
        topology=topology,
        sequence=sequence,
        sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        canonical_sha256=canonical_sequence_hash(
            sequence,
            circular=topology == "circular",
        ),
        source_input_id=f"sha256:{'b' * 64}",
        taxonomy={},
        metadata_extra={},
    )


def _truth_row(reference: Any, **changes: str) -> dict[str, str]:
    row = dict.fromkeys(validation_module.TRUTH_COLUMNS, "")
    row.update(
        {
            "sequence_id": f"frag-{'a' * 32}",
            "label": "virus",
            "source_accession_version": "ACCESSION.1",
            "source_genome_id": "genome-a",
            "source_content_group_id": f"sha256:{reference.digest}",
            "source_sequence_id": "sequence-a",
            "source_start": "0",
            "source_end": "4",
            "coordinate_system": "0-based-half-open",
            "strand": "+",
            "fragment_length": "4",
            "partition": "train",
            "view": "train",
            "release_date": "2026-01-01",
            "synthetic": "true",
        }
    )
    row.update(changes)
    return row


def test_validation_partition_reader_rejects_forged_truth_fields(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(validation_module, "validate_instance", lambda *_args, **_kwargs: None)
    reference = _reference(
        digest="d" * 64,
        accession_version="ACCESSION.1",
        release_date="2026-01-01",
        length_nt=4,
    )
    references = {reference.genome_id: reference}
    source = _fragment_source()
    sources = {source.sequence_id: source}
    fasta = tmp_path / "train.fasta.gz"
    truth = tmp_path / "train.truth.tsv.gz"
    base_contig = Contig(sequence_id=f"frag-{'a' * 32}", sequence="ACGT")

    def validate(
        row: dict[str, str],
        message: str,
        *,
        contig: Contig = base_contig,
        current_sources: dict[str, Any] = sources,
    ) -> None:
        write_fasta([contig], fasta, overwrite=True)
        write_tsv(truth, [row], validation_module.TRUTH_COLUMNS)
        with pytest.raises(IntegrityError, match=message):
            validation_module._read_partition(
                tmp_path,
                "train",
                sources=current_sources,
                references=references,
            )

    fasta.write_text("not FASTA\n", encoding="utf-8")
    write_tsv(truth, [_truth_row(reference)], validation_module.TRUTH_COLUMNS)
    with pytest.raises(IntegrityError, match="Cannot validate FASTA"):
        validation_module._read_partition(
            tmp_path,
            "train",
            sources=sources,
            references=references,
        )

    write_fasta([base_contig], fasta, overwrite=True)
    write_tsv(
        truth,
        [_truth_row(reference), _truth_row(reference)],
        validation_module.TRUTH_COLUMNS,
    )
    with pytest.raises(IntegrityError, match="duplicate truth sequence IDs"):
        validation_module._read_partition(
            tmp_path,
            "train",
            sources=sources,
            references=references,
        )

    cases = (
        (
            _truth_row(reference, sequence_id="readable-fragment-id"),
            "not an opaque CHIMERA ID",
            Contig(sequence_id="readable-fragment-id", sequence="ACGT"),
            sources,
        ),
        (_truth_row(reference, view="test"), "view must be 'train'", base_contig, sources),
        (
            _truth_row(reference, partition="test"),
            "semantic partition must be one of",
            base_contig,
            sources,
        ),
        (
            _truth_row(reference, label="other"),
            "label must be 'virus' or 'host'",
            base_contig,
            sources,
        ),
        (
            _truth_row(reference, fragment_length="3", source_end="3"),
            "fragment_length disagrees",
            base_contig,
            sources,
        ),
        (_truth_row(reference, strand="x"), "strand must be", base_contig, sources),
        (
            _truth_row(reference, source_sequence_id=""),
            "source identifiers must not be empty",
            base_contig,
            sources,
        ),
        (
            _truth_row(reference, source_accession_version="WRONG.1"),
            "accession disagrees",
            base_contig,
            sources,
        ),
        (
            _truth_row(reference, release_date="2025-01-01"),
            "release_date disagrees",
            base_contig,
            sources,
        ),
        (
            _truth_row(reference, source_content_group_id=f"sha256:{'0' * 64}"),
            "source content group digest is incorrect",
            base_contig,
            sources,
        ),
        (
            _truth_row(reference, source_end="5", fragment_length="5"),
            "fragment exceeds its source sequence length",
            Contig(sequence_id=f"frag-{'a' * 32}", sequence="ACGTA"),
            sources,
        ),
        (
            _truth_row(reference, source_start="1", source_end="5"),
            "linear coordinates exceed source length",
            base_contig,
            sources,
        ),
        (
            _truth_row(reference, coordinate_system="wrong"),
            "coordinate_system disagrees",
            base_contig,
            sources,
        ),
        (
            _truth_row(reference, synthetic="false"),
            "synthetic must be 'true'",
            base_contig,
            sources,
        ),
    )
    for row, message, contig, current_sources in cases:
        validate(row, message, contig=contig, current_sources=current_sources)

    circular_source = _fragment_source(topology="circular")
    validate(
        _truth_row(
            reference,
            source_start="4",
            source_end="8",
            coordinate_system="0-based-half-open-circular",
        ),
        "invalid circular source coordinates",
        current_sources={circular_source.sequence_id: circular_source},
    )
