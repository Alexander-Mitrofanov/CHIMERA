"""Tests for immutable sequence models and strict FASTA I/O."""

from __future__ import annotations

import gzip
import os
import stat
import tempfile
import unittest
from dataclasses import FrozenInstanceError, fields
from datetime import date, datetime
from pathlib import Path

from chimera.fasta import (
    DuplicateSequenceIdError,
    FastaDiscoveryError,
    FastaFormatError,
    discover_fasta_files,
    iter_fasta,
    read_fasta,
    write_fasta,
)
from chimera.models import (
    Contig,
    Fragment,
    Genome,
    GenomeMetadata,
    Label,
    canonical_sequence_hash,
    deterministic_genome_hash,
    deterministic_topology_agnostic_genome_hash,
    normalize_iupac_dna,
    reverse_complement,
    validate_stable_id,
)


class ModelTests(unittest.TestCase):
    def test_iupac_normalization_and_reverse_complement(self) -> None:
        self.assertEqual(
            normalize_iupac_dna(" acgt rysw\n kmbdhvn\t"),
            "ACGTRYSWKMBDHVN",
        )
        self.assertEqual(
            reverse_complement("ACGTRYSWKMBDHVN"),
            "NBDHVKMWSRYACGT",
        )

    def test_iupac_rejects_empty_rna_gaps_digits_and_unicode(self) -> None:
        for sequence, offending in (
            (" \n\t", "empty"),
            ("ACGU", "'U'"),
            ("AC-G", "'-'"),
            ("AC1G", "'1'"),
            ("ACßG", "'ß'"),
        ):
            with self.subTest(sequence=sequence), self.assertRaisesRegex(ValueError, offending):
                normalize_iupac_dna(sequence)
        with self.assertRaisesRegex(TypeError, "sequence must be a string"):
            normalize_iupac_dna(42)  # type: ignore[arg-type]

    def test_stable_id_contract(self) -> None:
        for identifier in ("NC_001.2", "genome-1", "x:y+z", "A" * 255):
            self.assertEqual(validate_stable_id(identifier), identifier)
        for identifier in ("", "-leading", "has space", "a/b", "a|b", "é", "A" * 256):
            with self.subTest(identifier=identifier), self.assertRaises(ValueError):
                validate_stable_id(identifier)
        with self.assertRaises(TypeError):
            validate_stable_id(42)  # type: ignore[arg-type]

    def test_contig_is_frozen_slotted_and_has_computed_properties(self) -> None:
        contig = Contig("seq-1", " acgn\n", "  example sequence  ", Path("in.fa"))
        self.assertEqual(contig.sequence, "ACGN")
        self.assertEqual(contig.description, "example sequence")
        self.assertEqual(contig.length, 4)
        self.assertEqual(contig.fasta_id, "seq-1")
        self.assertEqual(contig.digest, canonical_sequence_hash("ACGN"))
        self.assertFalse(hasattr(contig, "__dict__"))
        with self.assertRaises(FrozenInstanceError):
            contig.sequence = "AAAA"  # type: ignore[misc]

    def test_contig_rejects_multiline_description_and_non_path_source(self) -> None:
        with self.assertRaisesRegex(ValueError, "single line"):
            Contig("seq", "ACGT", "first\nsecond")
        with self.assertRaisesRegex(TypeError, "source_path"):
            Contig("seq", "ACGT", source_path="in.fa")  # type: ignore[arg-type]
        with self.assertRaisesRegex(TypeError, "description must be a string"):
            Contig("seq", "ACGT", description=1)  # type: ignore[arg-type]
        with self.assertRaisesRegex(TypeError, "release_date"):
            Contig("seq", "ACGT", release_date=datetime(2025, 1, 1))  # type: ignore[arg-type]
        with self.assertRaisesRegex(ValueError, "topology"):
            Contig("seq", "ACGT", topology="branched")  # type: ignore[arg-type]

    def test_canonical_sequence_hash_is_orientation_invariant(self) -> None:
        forward = "AAGCRYN"
        reverse = reverse_complement(forward)
        self.assertEqual(canonical_sequence_hash(forward), canonical_sequence_hash(reverse))
        self.assertNotEqual(canonical_sequence_hash(forward), canonical_sequence_hash("AAGCRYA"))

    def test_circular_hash_is_rotation_invariant_and_topology_is_domain_separated(self) -> None:
        sequence = "AACCGT"
        rotated = sequence[2:] + sequence[:2]
        circular = Contig("circular", sequence, topology="circular")
        rotated_circular = Contig("rotated", rotated, topology="circular")
        linear = Contig("linear", sequence, topology="linear")

        self.assertEqual(circular.digest, rotated_circular.digest)
        self.assertNotEqual(linear.digest, circular.digest)
        self.assertNotEqual(
            Genome("linear-genome", Label.VIRUS, (linear,)).digest,
            Genome("circular-genome", Label.VIRUS, (circular,)).digest,
        )
        with self.assertRaisesRegex(TypeError, "circular must be a boolean"):
            canonical_sequence_hash(sequence, circular=1)  # type: ignore[arg-type]

    def test_genome_hash_is_contig_order_name_and_orientation_invariant(self) -> None:
        first = Contig("a", "AACG")
        second = Contig("b", "TTAAN")
        renamed_reverse_first = Contig("renamed-1", reverse_complement(first.sequence))
        renamed_reverse_second = Contig("renamed-2", reverse_complement(second.sequence))

        original = Genome("g-1", Label.VIRUS, (first, second))
        equivalent = Genome(
            "g-2",
            "virus",
            (renamed_reverse_second, renamed_reverse_first),  # type: ignore[arg-type]
        )
        self.assertEqual(original.digest, equivalent.digest)
        self.assertEqual(original.length, 9)
        self.assertEqual(original.digest, deterministic_genome_hash(original.contigs))

        with_duplicate = Genome("g-3", Label.VIRUS, (first, second, Contig("duplicate", "AACG")))
        self.assertNotEqual(original.digest, with_duplicate.digest)

    def test_genome_digest_helpers_reject_invalid_contig_collections(self) -> None:
        for digest_function in (
            deterministic_genome_hash,
            deterministic_topology_agnostic_genome_hash,
        ):
            with (
                self.subTest(function=digest_function.__name__, value=None),
                self.assertRaisesRegex(TypeError, "iterable of Contig"),
            ):
                digest_function(None)  # type: ignore[arg-type]
            with (
                self.subTest(function=digest_function.__name__, value=()),
                self.assertRaisesRegex(ValueError, "at least one contig"),
            ):
                digest_function(())
            with (
                self.subTest(function=digest_function.__name__, value=("not-contig",)),
                self.assertRaisesRegex(TypeError, "only Contig"),
            ):
                digest_function(("not-contig",))  # type: ignore[arg-type]

    def test_genome_requires_nonempty_uniquely_named_contigs(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least one"):
            Genome("g", Label.HOST, ())
        with self.assertRaisesRegex(ValueError, "duplicates"):
            Genome(
                "g",
                Label.HOST,
                (Contig("same", "AAAA"), Contig("same", "CCCC")),
            )
        with self.assertRaisesRegex(ValueError, "virus.*host"):
            Genome("g", "bacterium", (Contig("c", "AAAA"),))  # type: ignore[arg-type]
        with self.assertRaisesRegex(TypeError, "iterable of Contig"):
            Genome("g", Label.HOST, "not-contigs")  # type: ignore[arg-type]
        with self.assertRaisesRegex(TypeError, "iterable of Contig"):
            Genome("g", Label.HOST, None)  # type: ignore[arg-type]
        with self.assertRaisesRegex(TypeError, "only Contig"):
            Genome("g", Label.HOST, ("not-a-contig",))  # type: ignore[arg-type]
        with self.assertRaisesRegex(TypeError, "GenomeMetadata"):
            Genome(
                "g",
                Label.HOST,
                (Contig("c", "AAAA"),),
                metadata={},  # type: ignore[arg-type]
            )

    def test_metadata_normalizes_taxonomy_and_uses_release_date(self) -> None:
        released = date(2025, 4, 3)
        metadata = GenomeMetadata(
            release_date=released,
            taxonomy=((" Family ", " Coronaviridae "), ("GENUS", "Betacoronavirus")),
            accession_version="NC_123456.2",
            extra=(("database", "NCBI"),),
        )
        self.assertEqual(metadata.taxon("family"), "Coronaviridae")
        self.assertEqual(metadata.taxon(" Genus "), "Betacoronavirus")
        self.assertIsNone(metadata.taxon("species"))
        self.assertEqual(metadata.release_date, released)
        self.assertEqual(metadata.deposited_at, released)
        self.assertIn("release_date", {item.name for item in fields(metadata)})
        self.assertNotIn("deposited_at", {item.name for item in fields(metadata)})
        self.assertFalse(hasattr(metadata, "__dict__"))

    def test_metadata_rejects_ambiguous_or_mutable_values(self) -> None:
        with self.assertRaisesRegex(ValueError, "duplicate key"):
            GenomeMetadata(taxonomy=(("Family", "A"), ("family", "B")))
        with self.assertRaisesRegex(TypeError, "release_date"):
            GenomeMetadata(release_date=datetime(2025, 1, 1))  # type: ignore[arg-type]
        with self.assertRaisesRegex(ValueError, "single line"):
            GenomeMetadata(extra=(("source", "one\ntwo"),))
        with self.assertRaisesRegex(TypeError, "iterable of key/value pairs"):
            GenomeMetadata(taxonomy="family=Alpha")  # type: ignore[arg-type]
        with self.assertRaisesRegex(TypeError, "iterable of key/value pairs"):
            GenomeMetadata(taxonomy=None)  # type: ignore[arg-type]
        with self.assertRaisesRegex(TypeError, "must be a .* pair"):
            GenomeMetadata(taxonomy=(("family",),))  # type: ignore[arg-type]
        with self.assertRaisesRegex(TypeError, "key must be a string"):
            GenomeMetadata(taxonomy=((1, "Alpha"),))  # type: ignore[arg-type]
        with self.assertRaisesRegex(ValueError, "key .* must start with a letter"):
            GenomeMetadata(taxonomy=(("1family", "Alpha"),))
        with self.assertRaisesRegex(ValueError, "must not be empty"):
            GenomeMetadata(taxonomy=(("family", "  "),))

    def test_fragment_validates_coordinates_and_keeps_label_out_of_fasta_id(self) -> None:
        fragment = Fragment(
            "frag-0001",
            " acgt ",
            "host",  # type: ignore[arg-type]
            "genome-1",
            "contig-1",
            10,
            14,
            "-",
            0,
        )
        self.assertEqual(fragment.label, Label.HOST)
        self.assertEqual(fragment.sequence, "ACGT")
        self.assertEqual(fragment.length, 4)
        self.assertEqual(fragment.fasta_id, "frag-0001")
        self.assertEqual(fragment.digest, canonical_sequence_hash("ACGT"))
        self.assertNotIn(fragment.label.value, fragment.fasta_id)
        self.assertFalse(hasattr(fragment, "__dict__"))

    def test_fragment_rejects_invalid_coordinates_strand_and_ordinal(self) -> None:
        base = {
            "fragment_id": "frag-1",
            "sequence": "ACGT",
            "label": Label.VIRUS,
            "genome_id": "g-1",
            "sequence_id": "c-1",
            "start": 0,
            "end": 4,
            "strand": "+",
            "ordinal": 0,
        }
        for replacement, message in (
            ({"start": -1, "end": 3}, "non-negative"),
            ({"end": 0}, "greater than start"),
            ({"end": 5}, "length"),
            ({"strand": "?"}, "strand"),
            ({"ordinal": -1}, "ordinal"),
            ({"start": True}, "integer"),
            ({"label": "bacterium"}, "virus.*host"),
        ):
            with self.subTest(replacement=replacement):
                values = base | replacement
                with self.assertRaisesRegex((TypeError, ValueError), message):
                    Fragment(**values)  # type: ignore[arg-type]


class FastaReaderTests(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self._temporary_directory.name)

    def tearDown(self) -> None:
        self._temporary_directory.cleanup()

    def test_discovery_is_recursive_sorted_case_insensitive_and_deduplicated(self) -> None:
        nested = self.root / "nested"
        nested.mkdir()
        first = self.root / "b.FASTA.GZ"
        with gzip.open(first, "wt", encoding="utf-8") as handle:
            handle.write(">b\nACGT\n")
        second = nested / "a.fna"
        second.write_text(">a\nTGCA\n", encoding="utf-8")
        (nested / "ignored.txt").write_text(">x\nAAAA\n", encoding="utf-8")

        discovered = discover_fasta_files((self.root, first, self.root))
        self.assertEqual(
            discovered,
            tuple(sorted((first.resolve(), second.resolve()), key=lambda p: p.as_posix())),
        )

    def test_plain_fasta_parses_multiline_iupac_description_and_source(self) -> None:
        source = self.root / "input.fa"
        source.write_text(
            "\n>seq-1 useful description\nacgt ry\nSWKM\n\n>seq-2\nNNNN\n",
            encoding="utf-8",
        )
        records = read_fasta(source)
        self.assertEqual([record.sequence_id for record in records], ["seq-1", "seq-2"])
        self.assertEqual(records[0].description, "useful description")
        self.assertEqual(records[0].sequence, "ACGTRYSWKM")
        self.assertEqual(records[0].source_path, source.resolve())
        self.assertEqual(records[1].sequence, "NNNN")

    def test_gzip_fasta_parses(self) -> None:
        source = self.root / "input.fasta.gz"
        with gzip.open(source, "wt", encoding="utf-8") as handle:
            handle.write(">compressed\nACGTN\n")
        [record] = read_fasta(source)
        self.assertEqual(record.sequence_id, "compressed")
        self.assertEqual(record.sequence, "ACGTN")

    def test_duplicate_ids_across_files_report_both_locations(self) -> None:
        first = self.root / "a.fa"
        second = self.root / "b.fa"
        first.write_text(">unique\nAAAA\n>duplicate\nCCCC\n", encoding="utf-8")
        second.write_text("\n>duplicate\nGGGG\n", encoding="utf-8")
        with self.assertRaises(DuplicateSequenceIdError) as caught:
            read_fasta(self.root)
        self.assertEqual(caught.exception.path, second.resolve())
        self.assertEqual(caught.exception.line_number, 2)
        self.assertIn("first seen at", str(caught.exception))
        self.assertIn(f"{first.resolve()}:3", str(caught.exception))

    def test_duplicate_ids_within_one_file_are_rejected(self) -> None:
        source = self.root / "duplicates.fa"
        source.write_text(">same\nAAAA\n>same\nCCCC\n", encoding="utf-8")
        with self.assertRaisesRegex(DuplicateSequenceIdError, "duplicate sequence ID"):
            tuple(iter_fasta(source))

    def test_format_errors_include_actionable_path_and_line(self) -> None:
        cases = (
            ("before.fa", "\nACGT\n", 2, "before the first FASTA header"),
            ("empty-header.fa", ">   \nACGT\n", 1, "header is empty"),
            ("empty-record.fa", ">empty\n>next\nAAAA\n", 1, "has no DNA sequence"),
            ("invalid-id.fa", ">has|pipe description\nAAAA\n", 1, "stable ID"),
            ("invalid-dna.fa", ">valid\nACGT\nAC-U\n", 3, "invalid sequence data"),
            ("empty.fa", "\n\n", 1, "contains no records"),
        )
        for filename, content, expected_line, message in cases:
            with self.subTest(filename=filename):
                source = self.root / filename
                source.write_text(content, encoding="utf-8")
                with self.assertRaises(FastaFormatError) as caught:
                    read_fasta(source)
                self.assertEqual(caught.exception.path, source.resolve())
                self.assertEqual(caught.exception.line_number, expected_line)
                self.assertIn(message, str(caught.exception))
                self.assertTrue(
                    str(caught.exception).startswith(f"{source.resolve()}:{expected_line}:")
                )

    def test_bad_gzip_and_invalid_utf8_have_line_aware_errors(self) -> None:
        bad_gzip = self.root / "bad.fa.gz"
        bad_gzip.write_bytes(b"not gzip")
        with self.assertRaises(FastaFormatError) as gzip_error:
            read_fasta(bad_gzip)
        self.assertEqual(gzip_error.exception.line_number, 1)
        self.assertIn("could not read", str(gzip_error.exception))

        invalid_utf8 = self.root / "bad.fa"
        invalid_utf8.write_bytes(b">valid\nACGT\n\xff\n")
        with self.assertRaises(FastaFormatError) as text_error:
            read_fasta(invalid_utf8)
        self.assertGreaterEqual(text_error.exception.line_number, 1)
        self.assertIn("UTF-8", str(text_error.exception))

    def test_empty_discovery_and_unsupported_explicit_file_are_actionable(self) -> None:
        self.assertEqual(discover_fasta_files(self.root), ())
        with self.assertRaisesRegex(FastaDiscoveryError, "no FASTA files found"):
            read_fasta(self.root)
        unsupported = self.root / "sequences.txt"
        unsupported.write_text(">x\nAAAA\n", encoding="utf-8")
        with self.assertRaisesRegex(FastaDiscoveryError, "unsupported FASTA filename"):
            discover_fasta_files(unsupported)
        with self.assertRaises(FileNotFoundError):
            discover_fasta_files(self.root / "missing.fa")


class FastaWriterTests(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self._temporary_directory.name)

    def tearDown(self) -> None:
        self._temporary_directory.cleanup()

    @staticmethod
    def fragment(*, fragment_id: str = "frag-0001", label: Label = Label.VIRUS) -> Fragment:
        return Fragment(
            fragment_id,
            "ACGTAC",
            label,
            "genome-1",
            "contig-1",
            2,
            8,
            "+",
            0,
        )

    def test_writer_wraps_contigs_and_round_trips(self) -> None:
        destination = self.root / "written.fasta"
        records = (
            Contig("seq-1", "ACGTACGT", "a description"),
            Contig("seq-2", "NNRY"),
        )
        returned = write_fasta(records, destination, line_width=3)
        self.assertEqual(returned, destination.resolve())
        self.assertEqual(
            destination.read_text(encoding="utf-8"),
            ">seq-1 a description\nACG\nTAC\nGT\n>seq-2\nNNR\nY\n",
        )
        round_tripped = read_fasta(destination)
        self.assertEqual(
            [(record.sequence_id, record.sequence) for record in round_tripped],
            [("seq-1", "ACGTACGT"), ("seq-2", "NNRY")],
        )

    def test_fragment_header_is_exactly_label_free_fragment_id(self) -> None:
        destination = self.root / "fragments.fa"
        fragments = (
            self.fragment(fragment_id="fragment-A", label=Label.VIRUS),
            self.fragment(fragment_id="fragment-B", label=Label.HOST),
        )
        write_fasta(fragments, destination)
        headers = [
            line
            for line in destination.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        ]
        self.assertEqual(headers, [">fragment-A", ">fragment-B"])
        self.assertNotIn("virus", "\n".join(headers).lower())
        self.assertNotIn("host", "\n".join(headers).lower())

    def test_gzip_output_is_readable_and_byte_reproducible(self) -> None:
        first = self.root / "first.fa.gz"
        second = self.root / "second.fa.gz"
        records = (Contig("seq", "ACGTN"),)
        write_fasta(records, first)
        write_fasta(records, second)
        self.assertEqual(first.read_bytes(), second.read_bytes())
        self.assertEqual(read_fasta(first)[0].sequence, "ACGTN")

    def test_writer_keeps_temp_private_and_publishes_shared_mode(self) -> None:
        destination = self.root / "permissions.fa"
        observed_temporary_mode: list[int] = []

        def records():
            temporary_files = list(self.root.glob(".permissions.fa.*.tmp"))
            self.assertEqual(len(temporary_files), 1)
            observed_temporary_mode.append(stat.S_IMODE(temporary_files[0].stat().st_mode))
            yield Contig("sequence", "ACGT")

        previous_umask = os.umask(0o077)
        try:
            write_fasta(records(), destination)
        finally:
            os.umask(previous_umask)

        self.assertEqual(observed_temporary_mode, [0o600])
        self.assertEqual(stat.S_IMODE(destination.stat().st_mode), 0o644)

    def test_existing_output_is_preserved_without_overwrite(self) -> None:
        destination = self.root / "existing.fa"
        destination.write_text("original\n", encoding="utf-8")
        with self.assertRaises(FileExistsError):
            write_fasta((Contig("new", "AAAA"),), destination)
        self.assertEqual(destination.read_text(encoding="utf-8"), "original\n")

        write_fasta((Contig("new", "AAAA"),), destination, overwrite=True)
        self.assertEqual(destination.read_text(encoding="utf-8"), ">new\nAAAA\n")

    def test_failed_stream_never_replaces_destination_or_leaves_temp_file(self) -> None:
        destination = self.root / "atomic.fa"
        destination.write_text("original\n", encoding="utf-8")

        def failing_records():
            yield Contig("first", "AAAA")
            raise RuntimeError("simulated generation failure")

        with self.assertRaisesRegex(RuntimeError, "simulated"):
            write_fasta(failing_records(), destination, overwrite=True)
        self.assertEqual(destination.read_text(encoding="utf-8"), "original\n")
        self.assertEqual(list(self.root.glob(".atomic.fa.*.tmp")), [])

    def test_writer_rejects_duplicate_ids_and_cleans_up(self) -> None:
        destination = self.root / "duplicates.fa"
        with self.assertRaisesRegex(ValueError, "duplicate FASTA output ID"):
            write_fasta((Contig("same", "AAAA"), Contig("same", "CCCC")), destination)
        self.assertFalse(destination.exists())
        self.assertEqual(list(self.root.glob(".duplicates.fa.*.tmp")), [])

    def test_writer_validates_destination_and_line_width(self) -> None:
        record = (Contig("seq", "AAAA"),)
        with self.assertRaises(FastaDiscoveryError):
            write_fasta(record, self.root / "output.txt")
        with self.assertRaises(ValueError):
            write_fasta(record, self.root / "output.fa", line_width=0)
        with self.assertRaises(TypeError):
            write_fasta(record, self.root / "output.fa", line_width=True)  # type: ignore[arg-type]
        with self.assertRaises(FileNotFoundError):
            write_fasta(record, self.root / "missing" / "output.fa")


if __name__ == "__main__":
    unittest.main()
