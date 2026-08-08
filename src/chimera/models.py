"""Immutable core sequence models and reproducible content digests.

Coordinates in this module are zero-based and half-open.  All biological
sequences are normalized to the uppercase, ungapped IUPAC DNA alphabet at
model boundaries so downstream algorithms never need to guess how to handle
case, whitespace, RNA, or alignment characters.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from datetime import date
from enum import StrEnum
from pathlib import Path
from typing import Final, Literal, cast

IUPAC_DNA_ALPHABET: Final[frozenset[str]] = frozenset("ACGTRYSWKMBDHVN")
"""The supported ungapped IUPAC DNA symbols."""

_IUPAC_DNA_INPUT_ALPHABET: Final[frozenset[str]] = frozenset(
    IUPAC_DNA_ALPHABET | {symbol.lower() for symbol in IUPAC_DNA_ALPHABET}
)
_COMPLEMENT: Final[dict[int, int]] = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_STABLE_ID_RE: Final[re.Pattern[str]] = re.compile(
    r"\A[A-Za-z0-9][A-Za-z0-9._:+-]{0,254}\Z", re.ASCII
)
_METADATA_KEY_RE: Final[re.Pattern[str]] = re.compile(r"\A[A-Za-z][A-Za-z0-9_.-]{0,63}\Z", re.ASCII)
Strand = Literal["+", "-"]
Topology = Literal["linear", "circular"]
CONTIG_DIGEST_ALGORITHM: Final[str] = "CHIMERA-CONTIG-SHA256-v2"
GENOME_DIGEST_ALGORITHM: Final[str] = "CHIMERA-GENOME-SHA256-v2"
TOPOLOGY_AGNOSTIC_GENOME_DIGEST_ALGORITHM: Final[str] = "CHIMERA-TOPOLOGY-AGNOSTIC-GENOME-SHA256-v1"


class Label(StrEnum):
    """Binary biological source label used throughout CHIMERA."""

    VIRUS = "virus"
    HOST = "host"


def validate_stable_id(value: str, *, field_name: str = "identifier") -> str:
    """Return *value* after validating the stable CHIMERA identifier syntax.

    Stable IDs are 1--255 ASCII characters, start with an alphanumeric
    character, and then contain only alphanumerics plus ``._:+-``.  This
    deliberately excludes whitespace and common serialization delimiters so
    an ID is represented identically in FASTA, TSV, JSON, and command lines.
    """

    if not isinstance(value, str):
        raise TypeError(f"{field_name} must be a string, got {type(value).__name__}")
    if not _STABLE_ID_RE.fullmatch(value):
        raise ValueError(
            f"{field_name} {value!r} is not a stable ID; use 1-255 ASCII "
            "characters, begin with a letter or digit, and use only letters, "
            "digits, '.', '_', ':', '+', or '-'"
        )
    return value


def normalize_iupac_dna(sequence: str) -> str:
    """Normalize DNA to uppercase IUPAC symbols and remove whitespace.

    ``U`` and gap symbols are rejected rather than silently reinterpreted.
    The error identifies the first offending character and its one-based
    position in the supplied string.
    """

    if not isinstance(sequence, str):
        raise TypeError(f"sequence must be a string, got {type(sequence).__name__}")

    compacted: list[str] = []
    for position, symbol in enumerate(sequence, start=1):
        if symbol.isspace():
            continue
        if symbol not in _IUPAC_DNA_INPUT_ALPHABET:
            raise ValueError(
                f"invalid IUPAC DNA symbol {symbol!r} at sequence position "
                f"{position}; allowed symbols are ACGTRYSWKMBDHVN"
            )
        compacted.append(symbol.upper())

    if not compacted:
        raise ValueError("DNA sequence is empty after removing whitespace")
    return "".join(compacted)


def reverse_complement(sequence: str) -> str:
    """Return the normalized IUPAC reverse complement of *sequence*."""

    normalized = normalize_iupac_dna(sequence)
    return normalized.translate(_COMPLEMENT)[::-1]


def _minimal_rotation(sequence: str) -> str:
    """Return the lexicographically minimal cyclic rotation in linear time."""

    doubled = sequence + sequence
    length = len(sequence)
    left, right, offset = 0, 1, 0
    while left < length and right < length and offset < length:
        left_symbol = doubled[left + offset]
        right_symbol = doubled[right + offset]
        if left_symbol == right_symbol:
            offset += 1
            continue
        if left_symbol > right_symbol:
            left = left + offset + 1
            if left == right:
                left += 1
        else:
            right = right + offset + 1
            if left == right:
                right += 1
        offset = 0
    start = min(left, right)
    return doubled[start : start + length]


def canonical_sequence_hash(sequence: str, *, circular: bool = False) -> str:
    """Return a versioned, strand- and topology-aware canonical SHA-256 digest.

    Linear records are invariant to strand orientation. Circular records are
    additionally invariant to rotation of the declared origin, using a
    linear-time minimal-rotation algorithm. The topology domain tag prevents a
    linear record and a circular record from being treated as the same source
    solely because their presented strings happen to match. Distinct ambiguity
    symbols remain distinct.
    """

    normalized = normalize_iupac_dna(sequence)
    reverse = normalized.translate(_COMPLEMENT)[::-1]
    if not isinstance(circular, bool):
        raise TypeError("circular must be a boolean")
    canonical = (
        min(_minimal_rotation(normalized), _minimal_rotation(reverse))
        if circular
        else min(normalized, reverse)
    )
    digest = hashlib.sha256()
    digest.update(f"{CONTIG_DIGEST_ALGORITHM}\0".encode("ascii"))
    digest.update(b"circular\0" if circular else b"linear\0")
    digest.update(len(canonical).to_bytes(8, byteorder="big", signed=False))
    digest.update(canonical.encode("ascii"))
    return digest.hexdigest()


def deterministic_genome_hash(contigs: tuple[Contig, ...]) -> str:
    """Return an order- and orientation-invariant SHA-256 genome digest.

    A genome is treated as a multiset of contig contents.  Duplicate contigs
    therefore remain significant, while contig names and input order do not.
    A domain marker and contig count make the serialization explicit and
    versionable.
    """

    try:
        immutable_contigs = tuple(contigs)
    except TypeError as error:
        raise TypeError("contigs must be an iterable of Contig objects") from error
    if not immutable_contigs:
        raise ValueError("a genome digest requires at least one contig")
    if not all(isinstance(contig, Contig) for contig in immutable_contigs):
        raise TypeError("contigs must contain only Contig objects")

    digests = sorted(bytes.fromhex(contig.digest) for contig in immutable_contigs)
    digest = hashlib.sha256()
    digest.update(f"{GENOME_DIGEST_ALGORITHM}\0".encode("ascii"))
    digest.update(len(digests).to_bytes(8, byteorder="big", signed=False))
    for contig_digest in digests:
        digest.update(contig_digest)
    return digest.hexdigest()


def deterministic_topology_agnostic_genome_hash(contigs: tuple[Contig, ...]) -> str:
    """Return a versioned exact raw/RC genome fingerprint for conflict audits.

    This deliberately ignores declared topology and does not normalize circular
    rotations. It catches identical presented sequence content (or its reverse
    complement) even when contradictory class labels carry different topology
    metadata. It is not the digest used for source grouping.
    """

    try:
        immutable_contigs = tuple(contigs)
    except TypeError as error:
        raise TypeError("contigs must be an iterable of Contig objects") from error
    if not immutable_contigs:
        raise ValueError("a genome digest requires at least one contig")
    if not all(isinstance(contig, Contig) for contig in immutable_contigs):
        raise TypeError("contigs must contain only Contig objects")

    digests = sorted(
        bytes.fromhex(canonical_sequence_hash(contig.sequence)) for contig in immutable_contigs
    )
    digest = hashlib.sha256()
    digest.update(f"{TOPOLOGY_AGNOSTIC_GENOME_DIGEST_ALGORITHM}\0".encode("ascii"))
    digest.update(len(digests).to_bytes(8, byteorder="big", signed=False))
    for contig_digest in digests:
        digest.update(contig_digest)
    return digest.hexdigest()


def _normalize_one_line_text(
    value: str,
    *,
    field_name: str,
    allow_empty: bool,
) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{field_name} must be a string, got {type(value).__name__}")
    if "\x00" in value or "\n" in value or "\r" in value:
        raise ValueError(f"{field_name} must be a single line without NUL bytes")
    normalized = value.strip()
    if not normalized and not allow_empty:
        raise ValueError(f"{field_name} must not be empty")
    return normalized


def _normalize_metadata_pairs(
    pairs: tuple[tuple[str, str], ...],
    *,
    field_name: str,
    lowercase_keys: bool,
) -> tuple[tuple[str, str], ...]:
    if isinstance(cast(object, pairs), (str, bytes)):
        raise TypeError(f"{field_name} must be an iterable of key/value pairs")
    try:
        immutable_pairs = tuple(pairs)
    except TypeError as error:
        raise TypeError(f"{field_name} must be an iterable of key/value pairs") from error

    normalized: list[tuple[str, str]] = []
    seen: set[str] = set()
    for index, pair in enumerate(immutable_pairs):
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            raise TypeError(f"{field_name}[{index}] must be a (key, value) pair")
        raw_key, raw_value = pair
        key = _normalize_one_line_text(
            raw_key, field_name=f"{field_name}[{index}] key", allow_empty=False
        )
        if lowercase_keys:
            key = key.lower()
        if not _METADATA_KEY_RE.fullmatch(key):
            raise ValueError(
                f"{field_name}[{index}] key {key!r} must start with a letter and "
                "contain only ASCII letters, digits, '.', '_', or '-'"
            )
        value = _normalize_one_line_text(
            raw_value,
            field_name=f"{field_name}[{index}] value",
            allow_empty=field_name == "extra",
        )
        if key in seen:
            raise ValueError(f"{field_name} contains duplicate key {key!r}")
        seen.add(key)
        normalized.append((key, value))
    return tuple(normalized)


@dataclass(frozen=True, slots=True)
class Contig:
    """A normalized FASTA record with a globally usable sequence ID."""

    sequence_id: str
    sequence: str
    description: str = ""
    source_path: Path | None = None
    accession_version: str | None = None
    release_date: date | None = None
    topology: Topology = "linear"
    taxonomy: tuple[tuple[str, str], ...] = ()
    metadata_extra: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "sequence_id",
            validate_stable_id(self.sequence_id, field_name="sequence_id"),
        )
        object.__setattr__(self, "sequence", normalize_iupac_dna(self.sequence))
        object.__setattr__(
            self,
            "description",
            _normalize_one_line_text(self.description, field_name="description", allow_empty=True),
        )
        if self.source_path is not None and not isinstance(self.source_path, Path):
            raise TypeError(
                f"source_path must be pathlib.Path or None, got {type(self.source_path).__name__}"
            )
        if self.accession_version is not None:
            object.__setattr__(
                self,
                "accession_version",
                validate_stable_id(self.accession_version, field_name="accession_version"),
            )
        if self.release_date is not None and type(self.release_date) is not date:
            raise TypeError("release_date must be datetime.date or None")
        if self.topology not in {"linear", "circular"}:
            raise ValueError("topology must be 'linear' or 'circular'")
        object.__setattr__(
            self,
            "taxonomy",
            _normalize_metadata_pairs(self.taxonomy, field_name="taxonomy", lowercase_keys=True),
        )
        object.__setattr__(
            self,
            "metadata_extra",
            _normalize_metadata_pairs(
                self.metadata_extra, field_name="extra", lowercase_keys=False
            ),
        )

    @property
    def length(self) -> int:
        """Normalized sequence length in nucleotides."""

        return len(self.sequence)

    @property
    def digest(self) -> str:
        """Canonical reverse-complement-aware content digest."""

        return canonical_sequence_hash(self.sequence, circular=self.topology == "circular")

    @property
    def fasta_id(self) -> str:
        """Identifier to place before the first whitespace in a FASTA header."""

        return self.sequence_id


@dataclass(frozen=True, slots=True)
class GenomeMetadata:
    """Immutable provenance and temporal/taxonomic metadata for one genome."""

    release_date: date | None = None
    taxonomy: tuple[tuple[str, str], ...] = ()
    accession_version: str | None = None
    extra: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        if self.release_date is not None and type(self.release_date) is not date:
            raise TypeError("release_date must be datetime.date or None")
        object.__setattr__(
            self,
            "taxonomy",
            _normalize_metadata_pairs(self.taxonomy, field_name="taxonomy", lowercase_keys=True),
        )
        if self.accession_version is not None:
            object.__setattr__(
                self,
                "accession_version",
                validate_stable_id(self.accession_version, field_name="accession_version"),
            )
        object.__setattr__(
            self,
            "extra",
            _normalize_metadata_pairs(self.extra, field_name="extra", lowercase_keys=False),
        )

    def taxon(self, rank: str) -> str | None:
        """Return the taxon name at *rank*, case-insensitively, if present."""

        normalized_rank = _normalize_one_line_text(
            rank, field_name="rank", allow_empty=False
        ).lower()
        for candidate_rank, taxon in self.taxonomy:
            if candidate_rank == normalized_rank:
                return taxon
        return None

    @property
    def deposited_at(self) -> date | None:
        """Read-only compatibility alias for :attr:`release_date`.

        New code and serialized metadata must use ``release_date``, meaning the
        NCBI first public release date.  This alias only supports migration of
        older split code.
        """

        return self.release_date


@dataclass(frozen=True, slots=True)
class Genome:
    """A labeled genome comprising one or more uniquely named contigs."""

    genome_id: str
    label: Label
    contigs: tuple[Contig, ...]
    metadata: GenomeMetadata = field(default_factory=GenomeMetadata)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "genome_id",
            validate_stable_id(self.genome_id, field_name="genome_id"),
        )
        try:
            object.__setattr__(self, "label", Label(self.label))
        except (TypeError, ValueError) as error:
            raise ValueError("label must be 'virus' or 'host'") from error
        if isinstance(cast(object, self.contigs), (str, bytes)):
            raise TypeError("contigs must be an iterable of Contig objects")
        try:
            immutable_contigs = tuple(self.contigs)
        except TypeError as error:
            raise TypeError("contigs must be an iterable of Contig objects") from error
        if not immutable_contigs:
            raise ValueError("genome must contain at least one contig")
        if not all(isinstance(contig, Contig) for contig in immutable_contigs):
            raise TypeError("contigs must contain only Contig objects")
        sequence_ids = [contig.sequence_id for contig in immutable_contigs]
        duplicate_ids = sorted(
            sequence_id for sequence_id in set(sequence_ids) if sequence_ids.count(sequence_id) > 1
        )
        if duplicate_ids:
            raise ValueError(
                "contig sequence IDs must be unique within a genome; duplicates: "
                + ", ".join(duplicate_ids)
            )
        object.__setattr__(self, "contigs", immutable_contigs)
        if not isinstance(self.metadata, GenomeMetadata):
            raise TypeError("metadata must be a GenomeMetadata object")

    @property
    def length(self) -> int:
        """Total genome length across all contigs."""

        return sum(contig.length for contig in self.contigs)

    @property
    def digest(self) -> str:
        """Deterministic content digest for the contig multiset."""

        return deterministic_genome_hash(self.contigs)


@dataclass(frozen=True, slots=True)
class Fragment:
    """A labeled fragment with auditable source coordinates.

    ``start`` is inclusive, ``end`` is exclusive, and ``ordinal`` is a
    zero-based stable position within the generating operation.  The FASTA ID
    is exactly ``fragment_id``; the class label is intentionally not encoded
    into it, preventing label leakage into downstream machine-learning inputs.
    """

    fragment_id: str
    sequence: str
    label: Label
    genome_id: str
    sequence_id: str
    start: int
    end: int
    strand: Strand
    ordinal: int

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "fragment_id",
            validate_stable_id(self.fragment_id, field_name="fragment_id"),
        )
        object.__setattr__(self, "sequence", normalize_iupac_dna(self.sequence))
        try:
            object.__setattr__(self, "label", Label(self.label))
        except (TypeError, ValueError) as error:
            raise ValueError("label must be 'virus' or 'host'") from error
        object.__setattr__(
            self,
            "genome_id",
            validate_stable_id(self.genome_id, field_name="genome_id"),
        )
        object.__setattr__(
            self,
            "sequence_id",
            validate_stable_id(self.sequence_id, field_name="sequence_id"),
        )
        for field_name in ("start", "end", "ordinal"):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"{field_name} must be an integer")
        if self.start < 0:
            raise ValueError("start must be non-negative")
        if self.end <= self.start:
            raise ValueError("end must be greater than start")
        if self.end - self.start != len(self.sequence):
            raise ValueError(
                "fragment sequence length must equal end - start for zero-based "
                "half-open coordinates"
            )
        if self.strand not in ("+", "-"):
            raise ValueError("strand must be '+' or '-'")
        if self.ordinal < 0:
            raise ValueError("ordinal must be non-negative")

    @property
    def length(self) -> int:
        """Fragment length in nucleotides."""

        return len(self.sequence)

    @property
    def digest(self) -> str:
        """Canonical reverse-complement-aware content digest."""

        return canonical_sequence_hash(self.sequence)

    @property
    def fasta_id(self) -> str:
        """Label-free identifier to write to a FASTA header."""

        return self.fragment_id


__all__ = [
    "CONTIG_DIGEST_ALGORITHM",
    "GENOME_DIGEST_ALGORITHM",
    "IUPAC_DNA_ALPHABET",
    "TOPOLOGY_AGNOSTIC_GENOME_DIGEST_ALGORITHM",
    "Contig",
    "Fragment",
    "Genome",
    "GenomeMetadata",
    "Label",
    "Strand",
    "canonical_sequence_hash",
    "deterministic_genome_hash",
    "deterministic_topology_agnostic_genome_hash",
    "normalize_iupac_dna",
    "reverse_complement",
    "validate_stable_id",
]
