"""Strict, streaming FASTA input and atomic FASTA output.

The reader accepts plain or gzip-compressed DNA FASTA files and reports the
source path and line number for malformed data.  IDs are required to be unique
across every file in one read operation.  The writer never appends biological
labels to fragment identifiers, which avoids a subtle source of label leakage
in machine-learning benchmarks.
"""

from __future__ import annotations

import errno
import gzip
import io
import os
import tempfile
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Final, TypeAlias

from .models import Contig, Fragment, normalize_iupac_dna, validate_stable_id

FASTA_SUFFIXES: Final[tuple[str, ...]] = (
    ".fa",
    ".fasta",
    ".fna",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
)
"""Recognized DNA FASTA filename suffixes, matched case-insensitively."""

Pathish: TypeAlias = str | os.PathLike[str]
FastaInputs: TypeAlias = Pathish | Iterable[Pathish]
FastaRecord: TypeAlias = Contig | Fragment

_PUBLISHED_FILE_MODE: Final = 0o644


class FastaError(ValueError):
    """Base class for FASTA discovery, format, and validation failures."""


class FastaDiscoveryError(FastaError):
    """Raised when input discovery cannot produce any readable FASTA files."""


class FastaFormatError(FastaError):
    """A malformed FASTA error tied to an exact source line."""

    def __init__(self, path: Path, line_number: int, message: str) -> None:
        self.path = path
        self.line_number = line_number
        self.message = message
        super().__init__(f"{path}:{line_number}: {message}")


class DuplicateSequenceIdError(FastaFormatError):
    """Raised when a sequence ID recurs within or across input FASTA files."""


def _is_fasta_path(path: Path) -> bool:
    filename = path.name.lower()
    return any(filename.endswith(suffix) for suffix in FASTA_SUFFIXES)


def _coerce_input_roots(inputs: FastaInputs) -> tuple[Path, ...]:
    roots: tuple[Path, ...]
    if isinstance(inputs, (str, os.PathLike)):
        roots = (Path(inputs),)
    else:
        try:
            roots = tuple(Path(item) for item in inputs)
        except TypeError as error:
            raise TypeError("FASTA inputs must be a path or an iterable of paths") from error
    if not roots:
        raise FastaDiscoveryError("no FASTA input paths were supplied")
    return roots


def discover_fasta_files(inputs: FastaInputs) -> tuple[Path, ...]:
    """Discover supported FASTA files beneath one or more paths.

    Directories are searched recursively.  Returned paths are resolved,
    de-duplicated (including repeated symlink targets), and sorted by their
    POSIX representation so file processing order is deterministic.  An empty
    directory is a valid discovery result; :func:`iter_fasta` turns that into
    an actionable error when records are actually requested.
    """

    discovered: dict[Path, Path] = {}
    for root in _coerce_input_roots(inputs):
        if not root.exists():
            raise FileNotFoundError(f"FASTA input does not exist: {root}")
        if root.is_file():
            if not _is_fasta_path(root):
                supported = ", ".join(FASTA_SUFFIXES)
                raise FastaDiscoveryError(
                    f"unsupported FASTA filename {root}; expected one of: {supported}"
                )
            resolved = root.resolve()
            discovered.setdefault(resolved, resolved)
            continue
        if not root.is_dir():
            raise FastaDiscoveryError(
                f"FASTA input is neither a regular file nor a directory: {root}"
            )
        for candidate in root.rglob("*"):
            if candidate.is_file() and _is_fasta_path(candidate):
                resolved = candidate.resolve()
                discovered.setdefault(resolved, resolved)

    return tuple(sorted(discovered.values(), key=lambda path: path.as_posix()))


def _parse_header(path: Path, line_number: int, line: str) -> tuple[str, str]:
    header = line[1:].strip()
    if not header:
        raise FastaFormatError(path, line_number, "FASTA header is empty after '>'")
    fields = header.split(maxsplit=1)
    sequence_id = fields[0]
    description = fields[1].strip() if len(fields) == 2 else ""
    try:
        validate_stable_id(sequence_id, field_name="sequence ID")
    except (TypeError, ValueError) as error:
        raise FastaFormatError(path, line_number, str(error)) from error
    if "\x00" in description:
        raise FastaFormatError(path, line_number, "FASTA description contains a NUL byte")
    return sequence_id, description


def _open_fasta_text(path: Path) -> io.TextIOBase:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8", newline=None)
    return path.open(mode="rt", encoding="utf-8", newline=None)


def _iter_fasta_file(
    path: Path,
    origins: dict[str, tuple[Path, int]],
) -> Iterator[Contig]:
    current_id: str | None = None
    current_description = ""
    current_header_line = 0
    sequence_chunks: list[str] = []
    line_number = 0
    record_count = 0

    def finish_record(*, next_header_line: int | None = None) -> Contig:
        nonlocal record_count
        assert current_id is not None
        if not sequence_chunks:
            detail = f"record {current_id!r} has no DNA sequence"
            if next_header_line is not None:
                detail += f" before the next header at line {next_header_line}"
            raise FastaFormatError(path, current_header_line, detail)
        try:
            contig = Contig(
                sequence_id=current_id,
                sequence="".join(sequence_chunks),
                description=current_description,
                source_path=path,
            )
        except (TypeError, ValueError) as error:
            raise FastaFormatError(path, current_header_line, str(error)) from error
        record_count += 1
        return contig

    try:
        with _open_fasta_text(path) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")
                if not line.strip():
                    continue
                if line.startswith(">"):
                    if current_id is not None:
                        yield finish_record(next_header_line=line_number)
                    sequence_id, description = _parse_header(path, line_number, line)
                    if sequence_id in origins:
                        first_path, first_line = origins[sequence_id]
                        raise DuplicateSequenceIdError(
                            path,
                            line_number,
                            f"duplicate sequence ID {sequence_id!r}; first seen at "
                            f"{first_path}:{first_line}",
                        )
                    origins[sequence_id] = (path, line_number)
                    current_id = sequence_id
                    current_description = description
                    current_header_line = line_number
                    sequence_chunks = []
                    continue

                if current_id is None:
                    raise FastaFormatError(
                        path,
                        line_number,
                        "DNA sequence appears before the first FASTA header; "
                        "headers must begin with '>' in column 1",
                    )
                try:
                    sequence_chunks.append(normalize_iupac_dna(line))
                except (TypeError, ValueError) as error:
                    raise FastaFormatError(
                        path,
                        line_number,
                        f"invalid sequence data in record {current_id!r}: {error}",
                    ) from error
    except FastaFormatError:
        raise
    except UnicodeError as error:
        raise FastaFormatError(
            path,
            max(line_number, 1),
            f"FASTA is not valid UTF-8 text: {error}",
        ) from error
    except OSError as error:
        raise FastaFormatError(
            path,
            max(line_number, 1),
            f"could not read FASTA data: {error}",
        ) from error

    if current_id is not None:
        yield finish_record()
    if record_count == 0:
        raise FastaFormatError(path, 1, "FASTA file contains no records")


def iter_fasta(inputs: FastaInputs) -> Iterator[Contig]:
    """Yield normalized contigs from all discovered inputs.

    Sequence IDs are checked globally for this complete iteration.  Callers
    that need all-or-nothing validation before processing should use
    :func:`read_fasta`, which materializes the result as a tuple.
    """

    files = discover_fasta_files(inputs)
    if not files:
        raise FastaDiscoveryError(
            "no FASTA files found; supported suffixes are " + ", ".join(FASTA_SUFFIXES)
        )
    origins: dict[str, tuple[Path, int]] = {}
    for path in files:
        yield from _iter_fasta_file(path, origins)


def read_fasta(inputs: FastaInputs) -> tuple[Contig, ...]:
    """Read and validate all discovered FASTA records into an immutable tuple."""

    return tuple(iter_fasta(inputs))


def _write_records(handle: io.TextIOBase, records: Iterable[FastaRecord], width: int) -> None:
    seen_ids: set[str] = set()
    for record_number, record in enumerate(records, start=1):
        if isinstance(record, Contig):
            fasta_id = record.sequence_id
            header = fasta_id
            if record.description:
                header += f" {record.description}"
        elif isinstance(record, Fragment):
            # Deliberately do not append label or provenance to fragment IDs.
            fasta_id = record.fragment_id
            header = fasta_id
        else:
            raise TypeError(
                f"FASTA record {record_number} must be Contig or Fragment, got "
                f"{type(record).__name__}"
            )
        if fasta_id in seen_ids:
            raise ValueError(f"duplicate FASTA output ID {fasta_id!r}")
        seen_ids.add(fasta_id)
        handle.write(f">{header}\n")
        for offset in range(0, len(record.sequence), width):
            handle.write(record.sequence[offset : offset + width])
            handle.write("\n")


def _write_temporary_fasta(
    descriptor: int,
    *,
    compressed: bool,
    records: Iterable[FastaRecord],
    width: int,
) -> None:
    with os.fdopen(descriptor, mode="wb", closefd=True) as raw_handle:
        if compressed:
            # An empty embedded filename and mtime=0 make gzip bytes reproducible.
            with (
                gzip.GzipFile(
                    filename="",
                    mode="wb",
                    fileobj=raw_handle,
                    mtime=0,
                ) as gzip_handle,
                io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="\n") as text_handle,
            ):
                _write_records(text_handle, records, width)
            raw_handle.flush()
            os.fchmod(raw_handle.fileno(), _PUBLISHED_FILE_MODE)
            os.fsync(raw_handle.fileno())
        else:
            with io.TextIOWrapper(raw_handle, encoding="utf-8", newline="\n") as text_handle:
                _write_records(text_handle, records, width)
                text_handle.flush()
                os.fchmod(raw_handle.fileno(), _PUBLISHED_FILE_MODE)
                os.fsync(raw_handle.fileno())


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        if error.errno in (errno.EINVAL, errno.ENOTSUP):
            return
        raise
    try:
        os.fsync(descriptor)
    except OSError as error:
        if error.errno not in (errno.EINVAL, errno.ENOTSUP):
            raise
    finally:
        os.close(descriptor)


def write_fasta(
    records: Iterable[FastaRecord],
    destination: Pathish,
    *,
    line_width: int = 80,
    overwrite: bool = False,
) -> Path:
    """Stream records to *destination* and commit them atomically.

    A temporary file is written and synchronized in the destination directory
    before a single atomic commit.  With ``overwrite=False`` (the default), an
    existing destination is never replaced, including if it appears during the
    write.  Files ending in ``.gz`` use deterministic gzip metadata.

    Fragment headers contain only ``fragment_id``.  Labels remain available in
    manifests/models but are never leaked through FASTA identifiers.
    """

    if isinstance(line_width, bool) or not isinstance(line_width, int):
        raise TypeError("line_width must be an integer")
    if line_width <= 0:
        raise ValueError("line_width must be greater than zero")
    destination_path = Path(destination)
    if not _is_fasta_path(destination_path):
        supported = ", ".join(FASTA_SUFFIXES)
        raise FastaDiscoveryError(
            f"unsupported FASTA output filename {destination_path}; expected one of: {supported}"
        )
    parent = destination_path.parent.resolve()
    if not parent.is_dir():
        raise FileNotFoundError(
            f"FASTA output directory does not exist or is not a directory: {parent}"
        )
    destination_path = parent / destination_path.name
    if destination_path.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing FASTA: {destination_path}")

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination_path.name}.", suffix=".tmp", dir=parent
    )
    temporary_path = Path(temporary_name)
    committed = False
    try:
        _write_temporary_fasta(
            descriptor,
            compressed=destination_path.name.lower().endswith(".gz"),
            records=records,
            width=line_width,
        )
        if overwrite:
            temporary_path.replace(destination_path)
        else:
            try:
                os.link(temporary_path, destination_path)
            except FileExistsError as error:
                raise FileExistsError(
                    f"refusing to overwrite existing FASTA: {destination_path}"
                ) from error
            temporary_path.unlink()
        committed = True
        _fsync_directory(parent)
        return destination_path
    finally:
        if not committed and temporary_path.exists():
            temporary_path.unlink()


__all__ = [
    "FASTA_SUFFIXES",
    "DuplicateSequenceIdError",
    "FastaDiscoveryError",
    "FastaError",
    "FastaFormatError",
    "discover_fasta_files",
    "iter_fasta",
    "read_fasta",
    "write_fasta",
]
