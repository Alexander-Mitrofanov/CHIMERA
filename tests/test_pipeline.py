from __future__ import annotations

import csv
import gzip
import json
import os
import random
import stat
import subprocess
from datetime import date
from enum import Enum
from pathlib import Path

import pytest

import chimera.pipeline as pipeline_module
from chimera.cli import main
from chimera.config import BenchmarkConfig, ConfigurationError, SimilarityBands, SplitKind
from chimera.errors import IntegrityError
from chimera.models import Contig, Fragment, Genome, Label
from chimera.pipeline import (
    BUNDLE_SCHEMA,
    _atomic_bundle_directory,
    _safe_target,
    generate_benchmark,
)
from chimera.references import ReferenceCatalog


def dna(seed: int, length: int = 320) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def publication_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    virus = tmp_path / "viruses.fna"
    host = tmp_path / "hosts.fna"
    metadata = tmp_path / "metadata.tsv"
    virus.write_text(
        f">v_old\n{dna(1)}\n>v_new\n{dna(2)}\n",
        encoding="utf-8",
    )
    host.write_text(
        f">h_old\n{dna(3)}\n>h_new\n{dna(4)}\n",
        encoding="utf-8",
    )
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\taccession_version\trelease_date\tfamily\n"
        "v_old\tv_old\tvirus\tNC_000001.1\t2018-01-01\tAlpha\n"
        "v_new\tv_new\tvirus\tNC_000002.1\t2022-01-01\tBeta\n"
        "h_old\th_old\thost\tGCF_000001.1\t2017-01-01\t\n"
        "h_new\th_new\thost\tGCF_000002.1\t2023-01-01\t\n",
        encoding="utf-8",
    )
    return virus, host, metadata


def suite_config(tmp_path: Path, *, output_name: str = "bundle") -> BenchmarkConfig:
    virus, host, metadata = publication_inputs(tmp_path)
    return BenchmarkConfig(
        virus_paths=(virus,),
        host_paths=(host,),
        metadata_path=metadata,
        output_dir=tmp_path / output_name,
        fragment_lengths=(31, 61),
        fragments_per_genome=8,
        similarity_k=15,
        sketch_size=200,
        seed=1234,
    )


def read_truth(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta_headers(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line[1:].strip() for line in handle if line.startswith(">")]


def test_one_command_generates_all_five_publication_protocols(tmp_path):
    config = suite_config(tmp_path)

    result = generate_benchmark(config)

    assert result.output_dir == config.output_dir
    assert (config.output_dir / "checksums.sha256").is_file()
    assert (config.output_dir / "manifest.json").is_file()
    assert (config.output_dir / "REPORT.md").is_file()
    assert (config.output_dir / "sequences.tsv").is_file()
    assert (config.output_dir / "source-sequences.fasta.gz").is_file()
    assert (config.output_dir / "schemas").is_dir()
    for kind in SplitKind:
        split_dir = config.output_dir / kind.directory_name
        assert (split_dir / "train.fasta.gz").is_file()
        assert (split_dir / "test.fasta.gz").is_file()
        assert (split_dir / "train.truth.tsv.gz").is_file()
        assert (split_dir / "test.truth.tsv.gz").is_file()
        manifest = json.loads((split_dir / "split.json").read_text())
        assert manifest["validation"]["status"] == "pass"
    similarity = config.output_dir / SplitKind.SIMILARITY.directory_name
    assert (similarity / "candidate_test.fasta.gz").is_file()
    for similarity_bin in (
        "high_similarity",
        "moderate_similarity",
        "low_similarity",
        "distant_detectable",
        "no_detectable_match",
    ):
        assert (similarity / "test_strata" / f"{similarity_bin}.fasta.gz").is_file()
        assert (similarity / "test_strata" / f"{similarity_bin}.truth.tsv.gz").is_file()


def test_generated_headers_do_not_leak_label_or_source(tmp_path):
    config = suite_config(tmp_path)
    generate_benchmark(config)
    split_dir = config.output_dir / SplitKind.GENOME.directory_name

    headers = read_fasta_headers(split_dir / "test.fasta.gz")
    truth = read_truth(split_dir / "test.truth.tsv.gz")

    assert headers
    assert set(headers) == {row["sequence_id"] for row in truth}
    assert all(header.startswith("frag-") for header in headers)
    assert all("virus" not in header and "host" not in header for header in headers)
    assert {row["label"] for row in truth} == {"virus", "host"}


def test_genome_temporal_and_taxonomy_invariants_are_materialized(tmp_path):
    config = suite_config(tmp_path)
    generate_benchmark(config)

    genome_rows = list(
        csv.DictReader(
            (config.output_dir / SplitKind.GENOME.directory_name / "assignments.tsv").open(),
            delimiter="\t",
        )
    )
    train = {row["genome_id"] for row in genome_rows if row["partition"] == "train"}
    test = {row["genome_id"] for row in genome_rows if row["partition"] == "test"}
    assert train.isdisjoint(test)

    temporal = json.loads(
        (config.output_dir / SplitKind.TEMPORAL.directory_name / "split.json").read_text()
    )
    assert temporal["parameters"]["temporal_semantics"] == "release-date-filtered retrospective"

    taxonomy_rows = list(
        csv.DictReader(
            (config.output_dir / SplitKind.TAXONOMY.directory_name / "assignments.tsv").open(),
            delimiter="\t",
        )
    )
    heldout = {row["taxon"] for row in taxonomy_rows if row["reason"] == "selected_taxon_holdout"}
    train_taxa = {row["taxon"] for row in taxonomy_rows if row["partition"] == "train"}
    assert heldout
    assert heldout.isdisjoint(train_taxa)


def test_force_rerun_is_byte_reproducible_except_execution_record(tmp_path):
    config = suite_config(tmp_path)
    generate_benchmark(config)
    first_checksums = (config.output_dir / "checksums.sha256").read_bytes()
    first_train = (
        config.output_dir / SplitKind.GENOME.directory_name / "train.fasta.gz"
    ).read_bytes()

    forced = BenchmarkConfig.from_mapping(
        {**config.as_manifest_dict(), "overwrite": True}, base_dir=tmp_path
    )
    generate_benchmark(forced)

    assert (config.output_dir / "checksums.sha256").read_bytes() == first_checksums
    assert (
        config.output_dir / SplitKind.GENOME.directory_name / "train.fasta.gz"
    ).read_bytes() == first_train


def test_force_refuses_to_delete_an_unrecognized_directory(tmp_path):
    config = suite_config(tmp_path)
    config.output_dir.mkdir()
    (config.output_dir / "user-data.txt").write_text("keep me", encoding="utf-8")
    forced = BenchmarkConfig.from_mapping(
        {**config.as_manifest_dict(), "overwrite": True}, base_dir=tmp_path
    )

    with pytest.raises(ConfigurationError, match="unrecognized"):
        generate_benchmark(forced)
    assert (config.output_dir / "user-data.txt").read_text() == "keep me"


def test_force_requires_the_exact_bundle_marker(tmp_path):
    config = suite_config(tmp_path)
    config.output_dir.mkdir()
    marker = config.output_dir / ".chimera-bundle"
    marker.write_text(f"{BUNDLE_SCHEMA}\ntrailing-data\n", encoding="utf-8")
    user_data = config.output_dir / "user-data.txt"
    user_data.write_text("keep me", encoding="utf-8")
    forced = BenchmarkConfig.from_mapping(
        {**config.as_manifest_dict(), "overwrite": True}, base_dir=tmp_path
    )

    with pytest.raises(ConfigurationError, match="must contain exactly"):
        generate_benchmark(forced)

    assert marker.read_text(encoding="utf-8").endswith("trailing-data\n")
    assert user_data.read_text(encoding="utf-8") == "keep me"


def test_force_refuses_a_marked_bundle_that_fails_validation(tmp_path):
    config = suite_config(tmp_path)
    generate_benchmark(config)
    report = config.output_dir / "REPORT.md"
    report.write_text(report.read_text(encoding="utf-8") + "\ntampered\n", encoding="utf-8")
    forced = BenchmarkConfig.from_mapping(
        {**config.as_manifest_dict(), "overwrite": True}, base_dir=tmp_path
    )

    with pytest.raises(ConfigurationError, match="invalid CHIMERA bundle"):
        generate_benchmark(forced)

    assert report.read_text(encoding="utf-8").endswith("\ntampered\n")


def test_safe_target_rejects_home_working_directory_and_their_ancestors(tmp_path):
    home = Path.home().resolve()
    working_directory = Path.cwd().resolve()
    protected = {Path("/"), home, working_directory, home.parent, working_directory.parent}

    for target in protected:
        with pytest.raises(ConfigurationError, match="broad/protected"):
            _safe_target(target)

    assert _safe_target(tmp_path / "bundle") == (tmp_path / "bundle").resolve()


def test_atomic_bundle_keeps_staging_private_and_publishes_shared_modes(tmp_path):
    target = tmp_path / "new-parent" / "bundle"
    previous_umask = os.umask(0o077)
    try:
        with _atomic_bundle_directory(target, overwrite=False) as staging:
            assert stat.S_IMODE(staging.stat().st_mode) == 0o700
            nested = staging / "nested"
            nested.mkdir()
            private_file = nested / "raw.txt"
            private_file.write_text("content", encoding="utf-8")
            assert stat.S_IMODE(nested.stat().st_mode) == 0o700
            assert stat.S_IMODE(private_file.stat().st_mode) == 0o600

        assert stat.S_IMODE(target.parent.stat().st_mode) == 0o755
        assert stat.S_IMODE(target.stat().st_mode) == 0o755
        assert stat.S_IMODE((target / "nested").stat().st_mode) == 0o755
        assert all(
            stat.S_IMODE(path.stat().st_mode) == 0o644
            for path in target.rglob("*")
            if path.is_file()
        )
    finally:
        os.umask(previous_umask)


def test_dry_run_resolves_every_plan_without_creating_output(tmp_path):
    config = suite_config(tmp_path)
    result = generate_benchmark(config, dry_run=True)
    assert result.dry_run is True
    assert not config.output_dir.exists()
    assert set(result.summary["resolved_plans"]) == {
        "genome",
        "similarity",
        "temporal",
        "taxonomy",
    }


def test_cli_suite_is_a_literal_one_command_workflow(tmp_path, capsys):
    virus, host, metadata = publication_inputs(tmp_path)
    output = tmp_path / "cli-bundle"
    status = main(
        [
            "suite",
            "--virus",
            str(virus),
            "--host",
            str(host),
            "--metadata",
            str(metadata),
            "--outdir",
            str(output),
            "--fragment-length",
            "31",
            "--fragments-per-genome",
            "6",
            "--similarity-k",
            "15",
            "--sketch-size",
            "200",
        ]
    )

    assert status == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["output_dir"] == str(output)
    assert (output / "manifest.json").is_file()


def test_source_sequence_copy_retains_per_contig_metadata_extra(tmp_path, monkeypatch):
    config = suite_config(tmp_path)
    metadata = config.metadata_path
    assert metadata is not None
    lines = metadata.read_text(encoding="utf-8").splitlines()
    metadata.write_text(
        "\n".join(
            [
                f"{lines[0]}\tsegment_note",
                *[f"{line}\tnote-{index}" for index, line in enumerate(lines[1:], 1)],
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    captured: dict[str, dict[str, str]] = {}
    real_write_fasta = pipeline_module.write_fasta

    def recording_write_fasta(records, destination, **kwargs):
        materialized = tuple(records)
        if Path(destination).name == "source-sequences.fasta.gz":
            captured.update(
                {contig.sequence_id: dict(contig.metadata_extra) for contig in materialized}
            )
        return real_write_fasta(materialized, destination, **kwargs)

    monkeypatch.setattr(pipeline_module, "write_fasta", recording_write_fasta)
    monkeypatch.setattr(pipeline_module, "validate_bundle", lambda _root: None)

    generate_benchmark(config)

    assert captured["v_old"]["segment_note"] == "note-1"
    assert captured["h_new"]["segment_note"] == "note-4"


def test_manifest_declares_linear_and_circular_coordinate_semantics(tmp_path, monkeypatch):
    config = suite_config(tmp_path)
    monkeypatch.setattr(pipeline_module, "validate_bundle", lambda _root: None)

    generate_benchmark(config)

    root_manifest = json.loads((config.output_dir / "manifest.json").read_text())
    assert root_manifest["data_model"]["coordinate_systems"] == {
        "linear": "0-based-half-open",
        "circular": "0-based-half-open-circular",
    }
    assert set(root_manifest["data_model"]["coordinate_semantics"]) == {
        "linear",
        "circular",
    }


def test_complete_staging_bundle_is_validated_before_atomic_publish(tmp_path, monkeypatch):
    config = suite_config(tmp_path)
    observed: list[tuple[Path, bool, bool]] = []

    def recording_validate_bundle(root):
        staging = Path(root)
        observed.append(
            (
                staging,
                config.output_dir.exists(),
                (staging / "checksums.sha256").is_file(),
            )
        )

    monkeypatch.setattr(pipeline_module, "validate_bundle", recording_validate_bundle)

    generate_benchmark(config)

    assert len(observed) == 1
    staging, output_existed, had_checksums = observed[0]
    assert staging.name.startswith(f".{config.output_dir.name}.staging-")
    assert output_existed is False
    assert had_checksums is True
    assert not staging.exists()
    assert config.output_dir.is_dir()


def test_jsonable_supports_public_provenance_types_and_rejects_unknown_values(tmp_path):
    class PlainEnum(Enum):
        VALUE = "plain-enum"

    value = {
        "date": date(2020, 1, 2),
        "path": tmp_path,
        "enum": SplitKind.GENOME,
        "bands": SimilarityBands(),
        "items": frozenset({1, 2}),
        "plain_enum": PlainEnum.VALUE,
    }
    serialized = pipeline_module._jsonable(value)
    assert serialized == {
        "bands": {"high": 0.9, "low": 0.3, "moderate": 0.7},
        "date": "2020-01-02",
        "enum": "genome",
        "items": [1, 2],
        "path": str(tmp_path),
        "plain_enum": "plain-enum",
    }
    with pytest.raises(TypeError, match="Cannot serialize provenance"):
        pipeline_module._jsonable(object())


def test_pipeline_directory_creation_rejects_file_parent(tmp_path):
    parent = tmp_path / "file"
    parent.write_text("content", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="Output parent is not a directory"):
        pipeline_module._ensure_public_directory(parent / "child")


def test_bundle_permission_publish_rejects_symlink(tmp_path):
    staging = tmp_path / "staging"
    staging.mkdir()
    target = staging / "target"
    target.write_text("content", encoding="utf-8")
    (staging / "link").symlink_to(target)
    with pytest.raises(IntegrityError, match="symbolic link"):
        pipeline_module._publish_bundle_permissions(staging)


def test_atomic_bundle_rejects_existing_target_without_force(tmp_path):
    target = tmp_path / "bundle"
    target.mkdir()
    with (
        pytest.raises(FileExistsError, match="use --force"),
        _atomic_bundle_directory(target, overwrite=False),
    ):
        pass


def test_atomic_bundle_rejects_preexisting_backup_path(tmp_path, monkeypatch):
    target = tmp_path / "bundle"
    target.mkdir()
    backup = tmp_path / f".bundle.backup-{os.getpid()}"
    backup.mkdir()
    monkeypatch.setattr(pipeline_module, "_require_valid_existing_bundle", lambda _path: None)

    with (
        pytest.raises(ConfigurationError, match="backup path already exists"),
        _atomic_bundle_directory(target, overwrite=True),
    ):
        pass
    assert target.is_dir()


def test_atomic_bundle_restores_target_if_staging_commit_fails(tmp_path, monkeypatch):
    target = tmp_path / "bundle"
    target.mkdir()
    original = target / "original.txt"
    original.write_text("keep", encoding="utf-8")
    monkeypatch.setattr(pipeline_module, "_require_valid_existing_bundle", lambda _path: None)
    real_replace = Path.replace

    def fail_staging_commit(path, destination):
        if path.name.startswith(".bundle.staging-") and Path(destination) == target:
            raise OSError("commit failed")
        return real_replace(path, destination)

    monkeypatch.setattr(Path, "replace", fail_staging_commit)
    with (
        pytest.raises(OSError, match="commit failed"),
        _atomic_bundle_directory(target, overwrite=True) as staging,
    ):
        (staging / "new.txt").write_text("new", encoding="utf-8")

    assert original.read_text(encoding="utf-8") == "keep"


def test_prepare_requires_two_genomes_per_label(tmp_path):
    virus = tmp_path / "virus.fna"
    host = tmp_path / "host.fna"
    virus.write_text(">v1\n" + dna(1, 80) + "\n", encoding="utf-8")
    host.write_text(">h1\n" + dna(2, 80) + "\n", encoding="utf-8")
    config = BenchmarkConfig(
        virus_paths=(virus,),
        host_paths=(host,),
        output_dir=tmp_path / "bundle",
        splits=(SplitKind.RANDOM,),
        fragment_lengths=(31,),
        fragments_per_genome=2,
    )
    with pytest.raises(ConfigurationError, match="At least two independent"):
        pipeline_module._prepare(config)


def test_prepare_rejects_fragment_length_unsupported_by_any_genome(tmp_path):
    virus_paths = []
    host_paths = []
    for label, destination in (("v", virus_paths), ("h", host_paths)):
        for index in (1, 2):
            path = tmp_path / f"{label}{index}.fna"
            path.write_text(f">{label}{index}\n{dna(index + (10 if label == 'h' else 0), 80)}\n")
            destination.append(path)
    config = BenchmarkConfig(
        virus_paths=tuple(virus_paths),
        host_paths=tuple(host_paths),
        output_dir=tmp_path / "bundle",
        splits=(SplitKind.RANDOM,),
        fragment_lengths=(100,),
        fragments_per_genome=2,
    )
    with pytest.raises(ConfigurationError, match="cannot emit 100-nt"):
        pipeline_module._prepare(config)


def _model_genome(genome_id: str, label: Label, sequence: str) -> Genome:
    return Genome(genome_id, label, (Contig(f"{genome_id}.1", sequence),))


def _model_fragment(
    fragment_id: str,
    genome: Genome,
    *,
    sequence: str = "ACGT",
) -> Fragment:
    return Fragment(
        fragment_id,
        sequence,
        genome.label,
        genome.genome_id,
        genome.contigs[0].sequence_id,
        0,
        len(sequence),
        "+",
        0,
    )


def test_split_validation_rejects_identifier_and_source_genome_leakage():
    genome = _model_genome("v", Label.VIRUS, "ACGT" * 20)
    fragment = _model_fragment("same", genome)
    with pytest.raises(IntegrityError, match="identifiers overlap"):
        pipeline_module._split_validation(
            SplitKind.GENOME, (fragment,), (fragment,), {"v": genome}, None
        )

    with pytest.raises(IntegrityError, match="source genomes leak"):
        pipeline_module._split_validation(
            SplitKind.GENOME,
            (_model_fragment("train", genome),),
            (_model_fragment("test", genome),),
            {"v": genome},
            None,
        )


def test_random_split_validation_requires_deliberate_source_overlap():
    train_genome = _model_genome("train", Label.VIRUS, "ACGT" * 20)
    test_genome = _model_genome("test", Label.VIRUS, "TGCA" * 20)
    with pytest.raises(IntegrityError, match="unexpectedly has no source-genome overlap"):
        pipeline_module._split_validation(
            SplitKind.RANDOM,
            (_model_fragment("train-fragment", train_genome),),
            (_model_fragment("test-fragment", test_genome),),
            {"train": train_genome, "test": test_genome},
            None,
        )


@pytest.mark.parametrize(("empty_side", "message"), [("train", "training"), ("test", "test")])
def test_split_validation_requires_each_label_on_each_side(empty_side, message):
    virus = _model_genome("virus", Label.VIRUS, dna(30, 80))
    host = _model_genome("host", Label.HOST, dna(31, 80))
    train = (_model_fragment("tv", virus), _model_fragment("th", host))
    test = (_model_fragment("xv", virus), _model_fragment("xh", host))
    if empty_side == "train":
        train = (_model_fragment("tv", virus),)
    else:
        test = (_model_fragment("xv", virus),)
    with pytest.raises(IntegrityError, match=rf"{message} partition lacks host"):
        pipeline_module._split_validation(
            SplitKind.RANDOM,
            train,
            test,
            {"virus": virus, "host": host},
            None,
        )


def test_split_validation_rejects_canonical_source_content_leakage():
    train_virus = _model_genome("tv", Label.VIRUS, dna(40, 80))
    test_virus = _model_genome("xv", Label.VIRUS, train_virus.contigs[0].sequence)
    train_host = _model_genome("th", Label.HOST, dna(41, 80))
    test_host = _model_genome("xh", Label.HOST, dna(42, 80))
    genomes = {g.genome_id: g for g in (train_virus, test_virus, train_host, test_host)}
    with pytest.raises(IntegrityError, match="genome content leaks"):
        pipeline_module._split_validation(
            SplitKind.GENOME,
            (_model_fragment("tvf", train_virus), _model_fragment("thf", train_host)),
            (_model_fragment("xvf", test_virus), _model_fragment("xhf", test_host)),
            genomes,
            None,
        )


def test_similarity_split_validation_rejects_exact_fragment_content_leakage():
    train_virus = _model_genome("tv", Label.VIRUS, dna(50, 80))
    test_virus = _model_genome("xv", Label.VIRUS, dna(51, 80))
    train_host = _model_genome("th", Label.HOST, dna(52, 80))
    test_host = _model_genome("xh", Label.HOST, dna(53, 80))
    genomes = {g.genome_id: g for g in (train_virus, test_virus, train_host, test_host)}
    with pytest.raises(IntegrityError, match="exact fragment content"):
        pipeline_module._split_validation(
            SplitKind.SIMILARITY,
            (
                _model_fragment("tvf", train_virus, sequence="AAAA"),
                _model_fragment("thf", train_host, sequence="CCCC"),
            ),
            (
                _model_fragment("xvf", test_virus, sequence="AAAA"),
                _model_fragment("xhf", test_host, sequence="GGGG"),
            ),
            genomes,
            None,
        )


def test_coordinate_overlap_counts_circular_origin_wrap():
    genome = Genome(
        "circular",
        Label.VIRUS,
        (Contig("circle.1", "ACGTACGTAA", topology="circular"),),
    )
    train = Fragment("train", "AAAC", Label.VIRUS, "circular", "circle.1", 8, 12, "+", 0)
    test = Fragment("test", "AAC", Label.VIRUS, "circular", "circle.1", 9, 12, "+", 1)
    assert (
        pipeline_module._test_fragments_with_coordinate_overlap(
            (train,), (test,), {genome.genome_id: genome}
        )
        == 1
    )


def test_report_requires_test_statistics():
    virus = _model_genome("v", Label.VIRUS, dna(60, 80))
    host = _model_genome("h", Label.HOST, dna(61, 80))
    prepared = pipeline_module._PreparedBenchmark(
        ReferenceCatalog((virus, host), ()),
        {},
    )
    config = BenchmarkConfig(
        virus_paths=(Path("v"),),
        host_paths=(Path("h"),),
        output_dir=Path("bundle"),
    )
    with pytest.raises(IntegrityError, match="lacks test statistics"):
        pipeline_module._report_markdown(config, prepared, {"genome": {}})


def test_generate_rejects_wrong_config_type_and_existing_output(tmp_path):
    with pytest.raises(TypeError, match="BenchmarkConfig"):
        generate_benchmark(object())

    config = suite_config(tmp_path)
    config.output_dir.mkdir()
    with pytest.raises(FileExistsError, match="use --force"):
        generate_benchmark(config)


def _external_similarity_config(tmp_path):
    base = suite_config(tmp_path, output_name="built-in-unused")
    prepared = pipeline_module._prepare(base)
    plan = prepared.plans[SplitKind.SIMILARITY]
    queries = [
        item.genome_id for item in plan.assignments if item.candidate_partition.value == "test"
    ]
    references = [
        item.genome_id for item in plan.assignments if item.candidate_partition.value == "train"
    ]
    table = tmp_path / "similarity.tsv"
    rows = [
        "query_genome_id\treference_genome_id\tsimilarity\tcoverage\tcoverage_definition\tmethod"
    ]
    rows.extend(
        f"{query}\t{reference}\t0.5\t1.0\taligned_fraction_shorter\tskani/0.3.1"
        for query in queries
        for reference in references
    )
    table.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return BenchmarkConfig.from_mapping(
        {
            **base.as_manifest_dict(),
            "output_dir": tmp_path / "external-bundle",
            "similarity_table": table,
        },
        base_dir=tmp_path,
    )


def test_external_similarity_table_is_snapshotted_and_content_addressed(tmp_path):
    config = _external_similarity_config(tmp_path)
    result = generate_benchmark(config)
    split_dir = result.output_dir / SplitKind.SIMILARITY.directory_name
    assert (
        split_dir / "external-similarity.tsv"
    ).read_bytes() == config.similarity_table.read_bytes()
    split_manifest = json.loads((split_dir / "split.json").read_text(encoding="utf-8"))
    source = split_manifest["parameters"]["similarity_source"]
    assert source.startswith("sha256:")
    assert split_manifest["parameters"]["similarity_table"] == source


def test_external_similarity_snapshot_read_failure_is_wrapped(tmp_path, monkeypatch):
    config = _external_similarity_config(tmp_path)
    table = config.similarity_table
    assert table is not None
    real_read_text = Path.read_text

    def fail_table(path, *args, **kwargs):
        if path == table:
            raise OSError("snapshot unavailable")
        return real_read_text(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", fail_table)
    with pytest.raises(ConfigurationError, match="Cannot snapshot external similarity table"):
        generate_benchmark(config)


@pytest.mark.parametrize("dirty", [False, True])
def test_git_provenance_accepts_valid_clean_and_dirty_repositories(tmp_path, monkeypatch, dirty):
    source_root = tmp_path / "source"
    (source_root / ".git").mkdir(parents=True)
    fake_file = source_root / "src/chimera/pipeline.py"
    fake_file.parent.mkdir(parents=True)
    monkeypatch.setattr(pipeline_module, "__file__", str(fake_file))
    revision = "a" * 40
    outputs = iter((revision + "\n", " M src/chimera/pipeline.py\n" if dirty else ""))

    def successful_run(*_args, **_kwargs):
        return subprocess.CompletedProcess([], 0, stdout=next(outputs), stderr="")

    monkeypatch.setattr(pipeline_module.subprocess, "run", successful_run)
    assert pipeline_module._git_provenance() == {
        "git_revision": revision,
        "git_dirty": dirty,
    }


def test_git_provenance_rejects_malformed_revision(tmp_path, monkeypatch):
    source_root = tmp_path / "source"
    (source_root / ".git").mkdir(parents=True)
    fake_file = source_root / "src/chimera/pipeline.py"
    fake_file.parent.mkdir(parents=True)
    monkeypatch.setattr(pipeline_module, "__file__", str(fake_file))
    outputs = iter(("not-a-revision\n", ""))
    monkeypatch.setattr(
        pipeline_module.subprocess,
        "run",
        lambda *_args, **_kwargs: subprocess.CompletedProcess(
            [], 0, stdout=next(outputs), stderr=""
        ),
    )
    assert pipeline_module._git_provenance() == {
        "git_revision": "unknown",
        "git_dirty": None,
    }


def test_git_provenance_handles_subprocess_failure(tmp_path, monkeypatch):
    source_root = tmp_path / "source"
    (source_root / ".git").mkdir(parents=True)
    fake_file = source_root / "src/chimera/pipeline.py"
    fake_file.parent.mkdir(parents=True)
    monkeypatch.setattr(pipeline_module, "__file__", str(fake_file))

    def fail(*_args, **_kwargs):
        raise subprocess.CalledProcessError(1, "git")

    monkeypatch.setattr(pipeline_module.subprocess, "run", fail)
    assert pipeline_module._git_provenance() == {
        "git_revision": "unknown",
        "git_dirty": None,
    }
