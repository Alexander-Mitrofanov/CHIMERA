from __future__ import annotations

import gzip
import os
import stat
from datetime import date
from pathlib import Path

import pytest

from chimera.errors import IntegrityError
from chimera.models import Contig, Fragment, Genome, GenomeMetadata, Label
from chimera.output import (
    REFERENCE_COLUMNS,
    TRUTH_COLUMNS,
    _atomic_binary_writer,
    _checksum_location,
    _ensure_public_directory,
    _resolved_bundle_root,
    _source_identifier,
    _validate_checksum_exclusions,
    _validate_relative_checksum_path,
    fragment_statistics,
    reference_rows,
    sequence_rows,
    summarize_split_truth,
    truth_rows,
    verify_checksums,
    write_checksums,
    write_json,
    write_text,
    write_tsv,
)


def source_genome() -> Genome:
    return Genome(
        "G1",
        Label.VIRUS,
        (Contig("NC_1.1", "ACGT" * 30),),
        GenomeMetadata(
            release_date=date(2020, 1, 2),
            accession_version="NC_1.1",
            taxonomy=(("family", "Exampleviridae"),),
        ),
    )


def fragments() -> tuple[Fragment, ...]:
    return (
        Fragment("f_opaque1", "ACGT", Label.VIRUS, "G1", "NC_1.1", 0, 4, "+", 0),
        Fragment("f_opaque2", "CGTA", Label.VIRUS, "G1", "NC_1.1", 1, 5, "-", 1),
    )


def test_json_and_text_writers_normalize_output(tmp_path):
    write_json(tmp_path / "value.json", {"z": 1, "a": 2})
    write_text(tmp_path / "value.txt", "a\r\nb")

    assert (tmp_path / "value.json").read_text() == '{\n  "a": 2,\n  "z": 1\n}\n'
    assert (tmp_path / "value.txt").read_bytes() == b"a\nb\n"


def test_atomic_writer_keeps_temp_private_and_publishes_shared_modes(tmp_path):
    destination = tmp_path / "new" / "nested" / "value.txt"
    previous_umask = os.umask(0o077)
    try:
        with _atomic_binary_writer(destination) as handle:
            temporary_files = list(destination.parent.glob(f".{destination.name}.*.tmp"))
            assert len(temporary_files) == 1
            assert stat.S_IMODE(temporary_files[0].stat().st_mode) == 0o600
            handle.write(b"publication data\n")
    finally:
        os.umask(previous_umask)

    assert destination.read_bytes() == b"publication data\n"
    assert stat.S_IMODE(destination.stat().st_mode) == 0o644
    assert stat.S_IMODE(destination.parent.stat().st_mode) == 0o755
    assert stat.S_IMODE(destination.parent.parent.stat().st_mode) == 0o755


def test_gzip_tsv_is_byte_reproducible_and_readable(tmp_path):
    rows = [{"a": "x", "b": 1}]
    first = tmp_path / "first.tsv.gz"
    second = tmp_path / "second.tsv.gz"
    write_tsv(first, rows, ("a", "b"))
    write_tsv(second, rows, ("a", "b"))

    assert first.read_bytes() == second.read_bytes()
    with gzip.open(first, "rt", encoding="utf-8") as handle:
        assert handle.read() == "a\tb\nx\t1\n"


def test_reference_and_truth_rows_preserve_provenance_but_ids_are_opaque():
    genome = source_genome()
    references = reference_rows([genome])
    truth = truth_rows(
        fragments(),
        partition="test",
        genomes={"G1": genome},
        assignment_details={
            "G1": {
                "similarity_bin": "moderate_similarity",
                "max_train_similarity": 0.8,
                "nearest_train_genome_id": "G0",
            }
        },
    )

    assert set(references[0]) == set(REFERENCE_COLUMNS)
    assert references[0]["release_date"] == "2020-01-02"
    assert set(truth[0]) == set(TRUTH_COLUMNS)
    assert truth[0]["sequence_id"] == "f_opaque1"
    assert "virus" not in str(truth[0]["sequence_id"])
    assert truth[0]["label"] == "virus"
    assert truth[0]["source_genome_id"] == "G1"
    assert truth[0]["source_content_group_id"] == f"sha256:{genome.digest}"
    assert truth[0]["coordinate_system"] == "0-based-half-open"
    assert truth[0]["partition"] == "test"
    assert truth[0]["view"] == "test"
    assert truth[0]["similarity_bin"] == "moderate_similarity"
    assert truth[0]["nearest_train_genome_id"] == "G0"


def test_truth_rows_distinguish_semantic_partition_from_serialized_view():
    genome = source_genome()

    truth = truth_rows(
        fragments(),
        partition={genome.genome_id: "excluded"},
        view="candidate_test",
        genomes={genome.genome_id: genome},
    )

    assert {row["partition"] for row in truth} == {"excluded"}
    assert {row["view"] for row in truth} == {"candidate_test"}


def test_truth_rows_use_exact_source_segment_release_date():
    first = Contig("segment-a", "ACGTACGT", release_date=date(2019, 1, 2))
    second = Contig("segment-b", "TGCATGCA", release_date=date(2020, 3, 4))
    genome = Genome(
        "segmented",
        Label.VIRUS,
        (first, second),
        GenomeMetadata(release_date=date(2020, 3, 4)),
    )
    segmented_fragments = (
        Fragment("frag-segment-a", "ACGT", Label.VIRUS, "segmented", "segment-a", 0, 4, "+", 0),
        Fragment("frag-segment-b", "TGCA", Label.VIRUS, "segmented", "segment-b", 0, 4, "+", 1),
    )

    rows = truth_rows(
        segmented_fragments,
        partition="test",
        genomes={genome.genome_id: genome},
    )

    dates_by_source = {row["source_sequence_id"]: row["release_date"] for row in rows}
    assert dates_by_source == {
        "segment-a": "2019-01-02",
        "segment-b": "2020-03-04",
    }


def test_fragment_statistics_report_effective_source_count():
    stats = fragment_statistics(fragments())
    assert stats["records"] == 2
    assert stats["bases"] == 8
    assert stats["source_genomes"] == 1
    assert stats["records_by_label"] == {"virus": 2}


def test_checksum_manifest_detects_corruption(tmp_path):
    write_text(tmp_path / "a.txt", "original")
    write_json(tmp_path / "b.json", {"ok": True})
    write_checksums(tmp_path)

    assert verify_checksums(tmp_path) == {"status": "pass", "files_checked": 2}

    write_text(tmp_path / "a.txt", "changed")
    with pytest.raises(IntegrityError, match=r"a\.txt"):
        verify_checksums(tmp_path)


def test_checksum_manifest_rejects_path_escape(tmp_path):
    write_text(tmp_path / "a.txt", "content")
    (tmp_path / "checksums.sha256").write_text("0" * 64 + "  ../outside\n")
    with pytest.raises(IntegrityError, match="escapes"):
        verify_checksums(tmp_path)


def test_checksum_manifest_must_be_nonempty(tmp_path):
    write_text(tmp_path / "a.txt", "content")
    write_text(tmp_path / "checksums.sha256", "")

    with pytest.raises(IntegrityError, match=r"manifest .* is empty"):
        verify_checksums(tmp_path)


def test_checksum_manifest_must_cover_exact_inventory(tmp_path):
    write_text(tmp_path / "a.txt", "a")
    write_text(tmp_path / "b.txt", "b")
    manifest = write_checksums(tmp_path)
    rows = manifest.read_text().splitlines()

    write_text(manifest, rows[0])
    with pytest.raises(IntegrityError, match=r"omitted file\(s\): b\.txt"):
        verify_checksums(tmp_path)

    write_text(manifest, "\n".join((*rows, f"{'0' * 64}  ghost.txt")))
    with pytest.raises(IntegrityError, match=r"unexpected path\(s\): ghost\.txt"):
        verify_checksums(tmp_path)


def test_checksum_manifest_rejects_duplicate_and_unsorted_rows(tmp_path):
    write_text(tmp_path / "a.txt", "a")
    write_text(tmp_path / "b.txt", "b")
    manifest = write_checksums(tmp_path)
    rows = manifest.read_text().splitlines()

    write_text(manifest, "\n".join((rows[0], rows[0], rows[1])))
    with pytest.raises(IntegrityError, match="duplicate path"):
        verify_checksums(tmp_path)

    write_text(manifest, "\n".join(reversed(rows)))
    with pytest.raises(IntegrityError, match="canonical sorted order"):
        verify_checksums(tmp_path)


@pytest.mark.parametrize(
    ("recorded_path", "message"),
    [
        ("./a.txt", "canonical POSIX"),
        ("nested/../a.txt", "canonical POSIX"),
        (r"nested\a.txt", "canonical POSIX"),
        ("/absolute/a.txt", "absolute checksum paths"),
        ("C:/absolute/a.txt", "absolute checksum paths"),
    ],
)
def test_checksum_manifest_rejects_noncanonical_paths(tmp_path, recorded_path, message):
    write_text(tmp_path / "a.txt", "content")
    write_text(tmp_path / "checksums.sha256", f"{'0' * 64}  {recorded_path}")

    with pytest.raises(IntegrityError, match=message):
        verify_checksums(tmp_path)


def test_checksum_writer_and_verifier_reject_symlinks_anywhere(tmp_path):
    target = write_text(tmp_path / "a.txt", "content")
    link = tmp_path / "linked.txt"
    link.symlink_to(target)

    with pytest.raises(IntegrityError, match="forbidden symbolic link"):
        write_checksums(tmp_path)

    link.unlink()
    write_checksums(tmp_path)
    (tmp_path / "execution.json").symlink_to(target)
    with pytest.raises(IntegrityError, match="forbidden symbolic link"):
        verify_checksums(tmp_path)


def test_checksum_inventory_excludes_only_manifest_and_execution_record(tmp_path):
    write_text(tmp_path / "a.txt", "content")
    execution = write_text(tmp_path / "execution.json", "volatile")
    manifest = write_checksums(tmp_path, destination=tmp_path / "integrity.sha256")

    assert manifest.read_text().split("  ", 1)[1].strip() == "a.txt"
    assert verify_checksums(tmp_path, manifest=manifest) == {
        "status": "pass",
        "files_checked": 1,
    }

    write_text(execution, "changed but deliberately unhashed")
    assert verify_checksums(tmp_path, manifest=manifest)["files_checked"] == 1

    with pytest.raises(ValueError, match=r"Only the checksum manifest and execution\.json"):
        write_checksums(tmp_path, destination=manifest, exclude=("a.txt",))


def test_checksum_writer_rejects_empty_bundle_inventory(tmp_path):
    with pytest.raises(IntegrityError, match="inventory is empty"):
        write_checksums(tmp_path)


def test_atomic_writer_removes_temporary_file_after_failure(tmp_path):
    destination = tmp_path / "failed.txt"
    with (
        pytest.raises(RuntimeError, match="stop"),
        _atomic_binary_writer(destination) as handle,
    ):
        handle.write(b"partial")
        raise RuntimeError("stop")

    assert not destination.exists()
    assert not list(tmp_path.glob(".failed.txt.*.tmp"))


def test_directory_creation_rejects_a_file_parent(tmp_path):
    parent = tmp_path / "not-a-directory"
    parent.write_text("content", encoding="utf-8")
    with pytest.raises(NotADirectoryError, match="Output parent"):
        _ensure_public_directory(parent / "child")


@pytest.mark.parametrize("digest", [None, "short", "G" * 64])
def test_source_identifier_requires_a_valid_available_digest(tmp_path, digest):
    source = tmp_path / "source.fna"
    mapping = {} if digest is None else {source.resolve(): digest}
    with pytest.raises(IntegrityError, match="No valid input digest"):
        _source_identifier(source, mapping)


def test_sequence_rows_keep_empty_source_identifier_for_library_contigs():
    rows = sequence_rows([source_genome()])
    assert rows[0]["source_input_id"] == ""


def test_truth_rows_reject_missing_mapped_partition_and_view():
    genome = source_genome()
    with pytest.raises(IntegrityError, match="No semantic partition"):
        truth_rows(
            fragments(),
            partition={},
            view="candidate_test",
            genomes={genome.genome_id: genome},
        )
    with pytest.raises(ValueError, match="view is required"):
        truth_rows(
            fragments(),
            partition={genome.genome_id: "test"},
            genomes={genome.genome_id: genome},
        )


@pytest.mark.parametrize("value", [True, "0.8"])
def test_truth_rows_reject_non_numeric_similarity(value):
    genome = source_genome()
    with pytest.raises(TypeError, match="max_train_similarity"):
        truth_rows(
            fragments(),
            partition="test",
            genomes={genome.genome_id: genome},
            assignment_details={genome.genome_id: {"max_train_similarity": value}},
        )


def test_resolved_bundle_root_rejects_symlink_missing_path_and_file(tmp_path):
    missing = tmp_path / "missing"
    with pytest.raises(IntegrityError, match="Cannot access bundle root"):
        _resolved_bundle_root(missing)

    regular = tmp_path / "file"
    regular.write_text("content", encoding="utf-8")
    with pytest.raises(IntegrityError, match="not a directory"):
        _resolved_bundle_root(regular)

    link = tmp_path / "link"
    link.symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(IntegrityError, match="symbolic link"):
        _resolved_bundle_root(link)


def test_checksum_location_accepts_relative_path_only_inside_root(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    path, relative = _checksum_location(tmp_path.resolve(), Path("nested/checksums.sha256"))
    assert path == (tmp_path / "nested/checksums.sha256").absolute()
    assert relative == "nested/checksums.sha256"

    monkeypatch.chdir(tmp_path.parent)
    with pytest.raises(IntegrityError, match="inside bundle root"):
        _checksum_location(tmp_path.resolve(), Path("outside.sha256"))


def test_checksum_path_and_exclusion_edge_contracts(tmp_path):
    with pytest.raises(IntegrityError, match="must not be empty"):
        _validate_relative_checksum_path("", path=tmp_path, line_number=4)
    with pytest.raises(ValueError, match="must not contain duplicates"):
        _validate_checksum_exclusions(("checksums.sha256", "checksums.sha256"), "checksums.sha256")


@pytest.mark.parametrize(
    ("manifest_text", "message"),
    [
        ("\n", "blank checksum row"),
        ("malformed\n", "malformed checksum row"),
        (f"{'G' * 64}  a.txt\n", "invalid SHA-256 digest"),
    ],
)
def test_checksum_manifest_rejects_malformed_rows(tmp_path, manifest_text, message):
    write_text(tmp_path / "a.txt", "content")
    (tmp_path / "checksums.sha256").write_text(manifest_text, encoding="utf-8")
    with pytest.raises(IntegrityError, match=message):
        verify_checksums(tmp_path)


def test_checksum_verifier_wraps_missing_manifest(tmp_path):
    write_text(tmp_path / "a.txt", "content")
    with pytest.raises(IntegrityError, match="Cannot read checksum manifest"):
        verify_checksums(tmp_path)


def test_checksum_writer_wraps_bundle_enumeration_failure(tmp_path, monkeypatch):
    write_text(tmp_path / "a.txt", "content")
    original_rglob = Path.rglob

    def failing_rglob(path, pattern):
        if path == tmp_path.resolve():
            raise OSError("enumeration failed")
        return original_rglob(path, pattern)

    monkeypatch.setattr(Path, "rglob", failing_rglob)
    with pytest.raises(IntegrityError, match="Cannot enumerate bundle"):
        write_checksums(tmp_path)


def test_split_truth_summary_is_sorted_and_counts_labels():
    assert summarize_split_truth(
        [
            {"partition": "test", "label": "virus"},
            {"partition": "train", "label": "host"},
            {"partition": "test", "label": "virus"},
        ]
    ) == {"test": {"virus": 2}, "train": {"host": 1}}
