from __future__ import annotations

from pathlib import Path

import pytest

import chimera.references as references_module
from chimera.config import DuplicatePolicy
from chimera.errors import InputError
from chimera.models import Contig, Label, reverse_complement
from chimera.references import (
    _load_labeled_contigs,
    _normalized_header,
    load_reference_catalog,
    read_metadata,
)


def write_fasta(path: Path, records: dict[str, str]) -> Path:
    path.write_text(
        "".join(f">{record_id}\n{sequence}\n" for record_id, sequence in records.items()),
        encoding="utf-8",
    )
    return path


def test_groups_segmented_genome_and_uses_latest_release_date(tmp_path):
    viruses = write_fasta(tmp_path / "virus.fna", {"segA": "ACGT" * 30, "segB": "TGCA" * 30})
    hosts = write_fasta(tmp_path / "host.fna", {"host1": "GATTACA" * 20})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\taccession_version\trelease_date\tfamily\n"
        "segA\tisolate1\tvirus\tNC_1.1\t2019-01-01\tFamilyA\n"
        "segB\tisolate1\tvirus\tNC_2.1\t2020-01-01\tFamilyA\n"
        "host1\thost1\thost\tGCF_1.1\t2018-01-01\t\n",
        encoding="utf-8",
    )

    catalog = load_reference_catalog((viruses,), (hosts,), metadata_path=metadata)
    grouped = catalog.by_id()["isolate1"]

    assert len(grouped.contigs) == 2
    assert grouped.metadata.release_date.isoformat() == "2020-01-01"
    assert grouped.metadata.taxon("family") == "FamilyA"
    assert grouped.metadata.accession_version is None


def test_grouped_genome_date_is_unknown_if_any_segment_date_is_missing(tmp_path):
    viruses = write_fasta(tmp_path / "virus.fna", {"segA": "ACGT" * 30, "segB": "TGCA" * 30})
    hosts = write_fasta(tmp_path / "host.fna", {"host1": "GATTACA" * 20})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\trelease_date\tfamily\n"
        "segA\tisolate1\tvirus\t2020-01-01\tFamilyA\n"
        "segB\tisolate1\tvirus\t\tFamilyA\n"
        "host1\thost1\thost\t2018-01-01\t\n",
        encoding="utf-8",
    )

    catalog = load_reference_catalog((viruses,), (hosts,), metadata_path=metadata)
    grouped = catalog.by_id()["isolate1"]

    assert grouped.metadata.release_date is None
    assert dict(grouped.metadata.extra)["missing_release_date_sequences"] == "segB"


def test_group_extra_requires_same_key_and_value_on_every_segment(tmp_path):
    viruses = write_fasta(tmp_path / "virus.fna", {"segA": "ACGT" * 30, "segB": "TGCA" * 30})
    hosts = write_fasta(tmp_path / "host.fna", {"host1": "GATTACA" * 20})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\tshared_note\tpartial_note\tconflicting_note\n"
        "segA\tisolate1\tvirus\tshared\tonly-segment-a\tvalue-a\n"
        "segB\tisolate1\tvirus\tshared\t\tvalue-b\n"
        "host1\thost1\thost\thost-value\t\thost-value\n",
        encoding="utf-8",
    )

    catalog = load_reference_catalog((viruses,), (hosts,), metadata_path=metadata)
    grouped = catalog.by_id()["isolate1"]
    group_extra = dict(grouped.metadata.extra)
    segment_extra = {contig.sequence_id: dict(contig.metadata_extra) for contig in grouped.contigs}

    assert group_extra["shared_note"] == "shared"
    assert "partial_note" not in group_extra
    assert "conflicting_note" not in group_extra
    assert segment_extra["segA"]["partial_note"] == "only-segment-a"
    assert segment_extra["segB"]["partial_note"] == ""
    assert segment_extra["segA"]["conflicting_note"] == "value-a"
    assert segment_extra["segB"]["conflicting_note"] == "value-b"


def test_blank_metadata_label_inherits_cli_source_class(tmp_path):
    viruses = write_fasta(tmp_path / "virus.fna", {"v1": "ACGT" * 30})
    hosts = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\nv1\tv1\t\nh1\th1\t\n",
        encoding="utf-8",
    )

    catalog = load_reference_catalog((viruses,), (hosts,), metadata_path=metadata)

    assert catalog.by_id()["v1"].label is Label.VIRUS
    assert catalog.by_id()["h1"].label is Label.HOST


def test_cross_class_content_collision_is_fatal_even_on_reverse_complement(tmp_path):
    sequence = "AACCGGTTACGA" * 10
    viruses = write_fasta(tmp_path / "virus.fna", {"v1": sequence})
    hosts = write_fasta(tmp_path / "host.fna", {"h1": reverse_complement(sequence)})

    with pytest.raises(InputError, match="Cross-class content conflict"):
        load_reference_catalog((viruses,), (hosts,))


def test_cross_class_exact_content_is_fatal_across_topology_declarations(tmp_path):
    sequence = "TACGAACCGGTT" * 10
    viruses = write_fasta(tmp_path / "virus.fna", {"v1": sequence})
    hosts = write_fasta(tmp_path / "host.fna", {"h1": sequence})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\ttopology\nv1\tv1\tvirus\tcircular\nh1\th1\thost\tlinear\n",
        encoding="utf-8",
    )

    with pytest.raises(InputError, match="topology-agnostic"):
        load_reference_catalog((viruses,), (hosts,), metadata_path=metadata)


def test_same_class_identical_sequence_with_different_topology_remains_distinct(tmp_path):
    sequence = "TACGAACCGGTT" * 10
    linear = write_fasta(tmp_path / "virus-linear.fna", {"v-linear": sequence})
    circular = write_fasta(tmp_path / "virus-circular.fna", {"v-circular": sequence})
    hosts = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\ttopology\n"
        "v-linear\tv-linear\tvirus\tlinear\n"
        "v-circular\tv-circular\tvirus\tcircular\n"
        "h1\th1\thost\tlinear\n",
        encoding="utf-8",
    )

    catalog = load_reference_catalog(
        (linear, circular),
        (hosts,),
        metadata_path=metadata,
        duplicate_policy=DuplicatePolicy.DROP,
    )

    retained = catalog.by_id()
    assert {"v-linear", "v-circular"} <= set(retained)
    assert retained["v-linear"].digest != retained["v-circular"].digest
    assert not catalog.exclusions


def test_same_class_duplicates_can_be_dropped_with_audit_row(tmp_path):
    virus_v2 = write_fasta(tmp_path / "virus-v2.fna", {"v2": "ACGT" * 30})
    virus_v1 = write_fasta(tmp_path / "virus-v1.fna", {"v1": "ACGT" * 30})
    hosts = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})

    catalog = load_reference_catalog(
        (virus_v2, virus_v1),
        (hosts,),
        duplicate_policy=DuplicatePolicy.DROP,
    )

    assert {item.genome_id for item in catalog.genomes} == {"v1", "h1"}
    assert catalog.exclusions[0].genome_id == "v2"
    assert catalog.exclusions[0].duplicate_of == "v1"


def test_same_class_duplicates_fail_by_default(tmp_path):
    virus_v1 = write_fasta(tmp_path / "virus-v1.fna", {"v1": "ACGT" * 30})
    virus_v2 = write_fasta(tmp_path / "virus-v2.fna", {"v2": "ACGT" * 30})
    hosts = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})
    with pytest.raises(InputError, match="Same-class duplicate"):
        load_reference_catalog((virus_v1, virus_v2), (hosts,))


def test_multi_record_fasta_requires_explicit_metadata(tmp_path):
    viruses = write_fasta(
        tmp_path / "virus.fna",
        {"segment-a": "ACGT" * 30, "segment-b": "TGCA" * 30},
    )
    hosts = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})

    with pytest.raises(InputError, match="Metadata is required for multi-record FASTA"):
        load_reference_catalog((viruses,), (hosts,))


def test_metadata_label_conflict_is_fatal(tmp_path):
    viruses = write_fasta(tmp_path / "virus.fna", {"v1": "ACGT" * 30})
    hosts = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("sequence_id\tlabel\nv1\thost\nh1\thost\n", encoding="utf-8")
    with pytest.raises(InputError, match="supplied via --virus"):
        load_reference_catalog((viruses,), (hosts,), metadata_path=metadata)


def test_missing_or_unused_metadata_rows_fail_without_silent_join(tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text("sequence_id\nunused\n", encoding="utf-8")
    viruses = write_fasta(tmp_path / "virus.fna", {"v1": "ACGT" * 30})
    hosts = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})
    with pytest.raises(InputError, match="has no row"):
        load_reference_catalog((viruses,), (hosts,), metadata_path=path)


def test_release_date_is_accepted_but_legacy_date_columns_are_rejected(tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text("sequence_id\trelease_date\nv1\t2020-01-02\n", encoding="utf-8")
    assert read_metadata(path)["v1"].release_date.isoformat() == "2020-01-02"

    for legacy_name in ("deposited_at", "create_date"):
        path.write_text(
            f"sequence_id\t{legacy_name}\nv1\t2020-01-02\n",
            encoding="utf-8",
        )
        with pytest.raises(InputError, match="ambiguous legacy date"):
            read_metadata(path)


def test_metadata_header_normalization_skips_none_and_rejects_duplicates():
    assert _normalized_header([None, " Sequence ID "]) == {"sequence_id": " Sequence ID "}
    with pytest.raises(InputError, match="duplicate column"):
        _normalized_header(["sequence-id", "Sequence ID"])


def test_metadata_open_failure_and_empty_table_are_actionable(tmp_path):
    with pytest.raises(InputError, match="Cannot read metadata table"):
        read_metadata(tmp_path / "missing.tsv")

    empty = tmp_path / "empty.tsv"
    empty.write_text("sequence_id\n", encoding="utf-8")
    with pytest.raises(InputError, match="has no data rows"):
        read_metadata(empty)


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ("genome_id\ng1\n", "missing required column"),
        ("sequence_id\tgenome_id\n\tg1\n", "sequence_id must not be empty"),
        ("sequence_id\nv1\nv1\n", "duplicate metadata sequence_id"),
        ("sequence_id\tlabel\nv1\tother\n", "label must be 'virus' or 'host'"),
        ("sequence_id\trelease_date\nv1\t2020/01/02\n", "release_date must use ISO"),
        ("sequence_id\ttopology\nv1\tbranched\n", "topology must be"),
    ],
)
def test_metadata_rejects_malformed_rows(tmp_path, payload, message):
    path = tmp_path / "metadata.tsv"
    path.write_text(payload, encoding="utf-8")
    with pytest.raises(InputError, match=message):
        read_metadata(path)


def test_metadata_csv_and_tax_prefix_are_supported(tmp_path):
    path = tmp_path / "metadata.csv"
    path.write_text(
        "sequence_id,genome_id,label,tax_family,custom\nv1,g1,VIRUS,Exampleviridae,receipt-1\n",
        encoding="utf-8",
    )
    record = read_metadata(path)["v1"]
    assert record.label is Label.VIRUS
    assert record.taxonomy == (("family", "Exampleviridae"),)
    assert record.extra == (("custom", "receipt-1"),)


def test_source_path_is_required_from_fasta_reader(monkeypatch):
    monkeypatch.setattr(
        references_module,
        "read_fasta",
        lambda _paths: (Contig("v1", "ACGT"),),
    )
    with pytest.raises(InputError, match="lacks a source path"):
        _load_labeled_contigs((Path("virus.fna"),), Label.VIRUS, None)


def test_grouped_genome_rejects_inconsistent_taxonomy(tmp_path):
    virus = write_fasta(tmp_path / "virus.fna", {"v1a": "ACGT" * 30, "v1b": "TGCA" * 30})
    host = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\tfamily\n"
        "v1a\tv1\tvirus\tAlpha\n"
        "v1b\tv1\tvirus\tBeta\n"
        "h1\th1\thost\t\n",
        encoding="utf-8",
    )
    with pytest.raises(InputError, match="inconsistent taxonomy"):
        load_reference_catalog((virus,), (host,), metadata_path=metadata)


def test_catalog_rejects_duplicate_sequence_ids_across_labels(tmp_path):
    virus = write_fasta(tmp_path / "virus.fna", {"same": "ACGT" * 30})
    host = write_fasta(tmp_path / "host.fna", {"same": "GATTACA" * 20})
    with pytest.raises(InputError, match="globally unique"):
        load_reference_catalog((virus,), (host,))


def test_catalog_rejects_unused_metadata_with_bounded_preview(tmp_path):
    virus = write_fasta(tmp_path / "virus.fna", {"v1": "ACGT" * 30})
    host = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})
    metadata = tmp_path / "metadata.tsv"
    extras = "".join(f"unused-{index}\n" for index in range(7))
    metadata.write_text(f"sequence_id\nv1\nh1\n{extras}", encoding="utf-8")
    with pytest.raises(InputError, match=r"7 sequence.*absent.*…"):
        load_reference_catalog((virus,), (host,), metadata_path=metadata)


def test_catalog_rejects_one_group_id_in_both_classes(tmp_path):
    virus = write_fasta(tmp_path / "virus.fna", {"v1": "ACGT" * 30})
    host = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\nv1\tshared\tvirus\nh1\tshared\thost\n",
        encoding="utf-8",
    )
    with pytest.raises(InputError, match="contradictory classes"):
        load_reference_catalog((virus,), (host,), metadata_path=metadata)


def test_invalid_group_identifier_is_wrapped_as_input_error(tmp_path):
    virus = write_fasta(tmp_path / "virus.fna", {"v1": "ACGT" * 30})
    host = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sequence_id\tgenome_id\tlabel\nv1\tbad id\tvirus\nh1\th1\thost\n",
        encoding="utf-8",
    )
    with pytest.raises(InputError, match="Invalid genome group"):
        load_reference_catalog((virus,), (host,), metadata_path=metadata)


def test_duplicate_drop_rejects_conflicting_taxonomy(tmp_path):
    first = write_fasta(tmp_path / "first.fna", {"v1": "ACGT" * 30})
    second = write_fasta(tmp_path / "second.fna", {"v2": "ACGT" * 30})
    host = write_fasta(tmp_path / "host.fna", {"h1": "GATTACA" * 20})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sequence_id\tlabel\tfamily\nv1\tvirus\tAlpha\nv2\tvirus\tBeta\nh1\thost\t\n",
        encoding="utf-8",
    )
    with pytest.raises(InputError, match="conflicting taxonomy"):
        load_reference_catalog(
            (first, second),
            (host,),
            metadata_path=metadata,
            duplicate_policy=DuplicatePolicy.DROP,
        )


def test_topology_aware_cross_class_digest_check_is_defensive(tmp_path, monkeypatch):
    virus = write_fasta(tmp_path / "virus.fna", {"v1": "ACGT" * 30})
    host = write_fasta(tmp_path / "host.fna", {"h1": "ACGT" * 30})
    fingerprints = iter(("a" * 64, "b" * 64))
    monkeypatch.setattr(
        references_module,
        "deterministic_topology_agnostic_genome_hash",
        lambda _contigs: next(fingerprints),
    )
    with pytest.raises(InputError, match="genome digest"):
        load_reference_catalog((virus,), (host,))
