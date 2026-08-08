"""Deterministic, atomic writers for CHIMERA benchmark bundles."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import os
import tempfile
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any

from .errors import IntegrityError
from .models import Fragment, Genome
from .similarity import format_similarity_value

_PRIVATE_DIRECTORY_MODE = 0o700
_PUBLISHED_DIRECTORY_MODE = 0o755
_PUBLISHED_FILE_MODE = 0o644


def _ensure_public_directory(path: Path) -> None:
    """Create missing path components privately, then publish each as ``0755``."""

    missing: list[Path] = []
    cursor = path
    while not cursor.exists():
        missing.append(cursor)
        parent = cursor.parent
        if parent == cursor:
            break
        cursor = parent
    if cursor.exists() and not cursor.is_dir():
        raise NotADirectoryError(f"Output parent is not a directory: {cursor}")
    for directory in reversed(missing):
        try:
            directory.mkdir(mode=_PRIVATE_DIRECTORY_MODE)
        except FileExistsError as error:
            if not directory.is_dir():
                raise NotADirectoryError(
                    f"Output parent is not a directory: {directory}"
                ) from error
        else:
            directory.chmod(_PUBLISHED_DIRECTORY_MODE)


def sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Return a lowercase SHA-256 digest without loading a file into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


@contextmanager
def _atomic_binary_writer(path: Path) -> Iterator[io.BufferedWriter]:
    _ensure_public_directory(path.parent)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    committed = False
    try:
        with os.fdopen(descriptor, "wb") as handle:
            yield handle
            handle.flush()
            os.fchmod(handle.fileno(), _PUBLISHED_FILE_MODE)
            os.fsync(handle.fileno())
        temporary.replace(path)
        committed = True
    finally:
        if not committed and temporary.exists():
            temporary.unlink()


def write_json(path: Path, value: Any) -> Path:
    """Write canonical, human-readable JSON atomically."""

    payload = json.dumps(
        value,
        ensure_ascii=False,
        indent=2,
        sort_keys=True,
        allow_nan=False,
    )
    with _atomic_binary_writer(path) as handle:
        handle.write(payload.encode("utf-8"))
        handle.write(b"\n")
    return path


def write_text(path: Path, value: str) -> Path:
    """Write normalized UTF-8 text atomically."""

    normalized = value.replace("\r\n", "\n").replace("\r", "\n")
    if normalized and not normalized.endswith("\n"):
        normalized += "\n"
    with _atomic_binary_writer(path) as handle:
        handle.write(normalized.encode("utf-8"))
    return path


def write_tsv(
    path: Path,
    rows: Iterable[Mapping[str, object]],
    columns: Sequence[str],
    *,
    compressed: bool | None = None,
) -> Path:
    """Write a stable RFC-4180-style tab-separated table atomically.

    Gzip output is reproducible (mtime zero and no embedded filename).
    """

    use_gzip = path.name.endswith(".gz") if compressed is None else compressed
    with _atomic_binary_writer(path) as raw:
        if use_gzip:
            binary: Any = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
        else:
            binary = raw
        with binary if use_gzip else _null_context(binary):
            text = io.TextIOWrapper(binary, encoding="utf-8", newline="", write_through=True)
            writer = csv.DictWriter(
                text,
                fieldnames=list(columns),
                delimiter="\t",
                lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({column: row.get(column, "") for column in columns})
            text.flush()
            if use_gzip:
                text.detach()
    return path


@contextmanager
def _null_context(value: Any) -> Iterator[Any]:
    yield value


def _source_identifier(
    path: Path,
    source_digests: Mapping[Path, str] | None,
) -> str:
    resolved = path.resolve()
    if source_digests is None:
        # Library-only row construction may not have the input bytes available.
        # Hashing the locator avoids disclosing it; publication generation always
        # supplies content digests and therefore emits the stronger ``sha256:`` form.
        locator = hashlib.sha256(str(resolved).encode("utf-8")).hexdigest()
        return f"locator-sha256:{locator}"
    digest = source_digests.get(resolved)
    if (
        digest is None
        or len(digest) != 64
        or any(char not in "0123456789abcdef" for char in digest)
    ):
        raise IntegrityError(f"No valid input digest is available for source file {resolved}")
    return f"sha256:{digest}"


def reference_rows(
    genomes: Iterable[Genome],
    *,
    source_digests: Mapping[Path, str] | None = None,
) -> list[dict[str, object]]:
    """Build the versioned reference inventory rows."""

    rows: list[dict[str, object]] = []
    for genome in sorted(genomes, key=lambda item: item.genome_id):
        source_paths = sorted(
            {
                _source_identifier(contig.source_path, source_digests)
                for contig in genome.contigs
                if contig.source_path is not None
            }
        )
        rows.append(
            {
                "genome_id": genome.genome_id,
                "label": genome.label.value,
                "accession_version": genome.metadata.accession_version or "",
                "release_date": (
                    genome.metadata.release_date.isoformat()
                    if genome.metadata.release_date is not None
                    else ""
                ),
                "sequence_ids": json.dumps(
                    [contig.sequence_id for contig in genome.contigs],
                    ensure_ascii=False,
                    separators=(",", ":"),
                ),
                "contig_count": len(genome.contigs),
                "length_nt": genome.length,
                "sha256": genome.digest,
                "source_input_ids": json.dumps(
                    source_paths,
                    ensure_ascii=False,
                    separators=(",", ":"),
                ),
                "taxonomy": json.dumps(
                    dict(genome.metadata.taxonomy),
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                "metadata_extra": json.dumps(
                    dict(genome.metadata.extra),
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ),
            }
        )
    return rows


REFERENCE_COLUMNS = (
    "genome_id",
    "label",
    "accession_version",
    "release_date",
    "sequence_ids",
    "contig_count",
    "length_nt",
    "sha256",
    "source_input_ids",
    "taxonomy",
    "metadata_extra",
)

SEQUENCE_COLUMNS = (
    "sequence_id",
    "genome_id",
    "label",
    "accession_version",
    "release_date",
    "topology",
    "length_nt",
    "sha256",
    "canonical_sha256",
    "source_input_id",
    "taxonomy",
    "metadata_extra",
)


def sequence_rows(
    genomes: Iterable[Genome],
    *,
    source_digests: Mapping[Path, str] | None = None,
) -> list[dict[str, object]]:
    """Build the exact per-contig inventory used to authenticate fragment truth."""

    rows: list[dict[str, object]] = []
    for genome in sorted(genomes, key=lambda item: item.genome_id):
        rows.extend(
            (
                {
                    "sequence_id": contig.sequence_id,
                    "genome_id": genome.genome_id,
                    "label": genome.label.value,
                    "accession_version": contig.accession_version or "",
                    "release_date": contig.release_date.isoformat() if contig.release_date else "",
                    "topology": contig.topology,
                    "length_nt": contig.length,
                    "sha256": hashlib.sha256(contig.sequence.encode("ascii")).hexdigest(),
                    "canonical_sha256": contig.digest,
                    "source_input_id": (
                        _source_identifier(contig.source_path, source_digests)
                        if contig.source_path
                        else ""
                    ),
                    "taxonomy": json.dumps(
                        dict(contig.taxonomy),
                        ensure_ascii=False,
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                    "metadata_extra": json.dumps(
                        dict(contig.metadata_extra),
                        ensure_ascii=False,
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                }
            )
            for contig in sorted(genome.contigs, key=lambda item: item.sequence_id)
        )
    return rows


TRUTH_COLUMNS = (
    "sequence_id",
    "label",
    "source_accession_version",
    "source_genome_id",
    "source_content_group_id",
    "source_sequence_id",
    "source_start",
    "source_end",
    "coordinate_system",
    "strand",
    "fragment_length",
    "partition",
    "view",
    "similarity_bin",
    "max_train_similarity",
    "nearest_train_genome_id",
    "release_date",
    "synthetic",
)


def truth_rows(
    fragments: Iterable[Fragment],
    *,
    partition: str | Mapping[str, str],
    view: str | None = None,
    genomes: Mapping[str, Genome],
    assignment_details: Mapping[str, Mapping[str, object]] | None = None,
) -> list[dict[str, object]]:
    """Build one unambiguous truth row per emitted fragment."""

    details = assignment_details or {}
    contigs = {
        (genome.genome_id, contig.sequence_id): contig
        for genome in genomes.values()
        for contig in genome.contigs
    }
    rows: list[dict[str, object]] = []
    for fragment in fragments:
        genome = genomes[fragment.genome_id]
        source_contig = contigs[(fragment.genome_id, fragment.sequence_id)]
        assignment = details.get(fragment.genome_id, {})
        if isinstance(partition, str):
            semantic_partition = partition
            view_name = view or partition
        else:
            try:
                semantic_partition = partition[fragment.genome_id]
            except KeyError as error:
                raise IntegrityError(
                    f"No semantic partition was supplied for genome {fragment.genome_id!r}"
                ) from error
            if view is None:
                raise ValueError("view is required when partition is supplied as a mapping")
            view_name = view
        similarity = assignment.get("max_train_similarity")
        if similarity is None:
            similarity_text = ""
        elif isinstance(similarity, (int, float)) and not isinstance(similarity, bool):
            similarity_text = format_similarity_value(float(similarity))
        else:
            raise TypeError("max_train_similarity must be a number or None")
        rows.append(
            {
                "sequence_id": fragment.fragment_id,
                "label": fragment.label.value,
                "source_accession_version": source_contig.accession_version or "",
                "source_genome_id": fragment.genome_id,
                "source_content_group_id": f"sha256:{genome.digest}",
                "source_sequence_id": fragment.sequence_id,
                "source_start": fragment.start,
                "source_end": fragment.end,
                "coordinate_system": (
                    "0-based-half-open-circular"
                    if source_contig.topology == "circular"
                    else "0-based-half-open"
                ),
                "strand": fragment.strand,
                "fragment_length": fragment.length,
                "partition": semantic_partition,
                "view": view_name,
                "similarity_bin": assignment.get("similarity_bin") or "",
                "max_train_similarity": similarity_text,
                "nearest_train_genome_id": assignment.get("nearest_train_genome_id") or "",
                "release_date": (
                    source_contig.release_date.isoformat()
                    if source_contig.release_date is not None
                    else ""
                ),
                "synthetic": "true",
            }
        )
    return rows


def fragment_statistics(fragments: Iterable[Fragment]) -> dict[str, object]:
    """Compute transparent composition summaries for a fragment collection."""

    materialized = tuple(fragments)
    by_label = Counter(fragment.label.value for fragment in materialized)
    by_length = Counter(str(fragment.length) for fragment in materialized)
    by_label_and_length = Counter(
        f"{fragment.label.value}:{fragment.length}" for fragment in materialized
    )
    by_genome = Counter(fragment.genome_id for fragment in materialized)
    total_bases = sum(fragment.length for fragment in materialized)
    gc_bases = sum(
        fragment.sequence.count("G") + fragment.sequence.count("C") for fragment in materialized
    )
    ambiguous_bases = sum(
        sum(base not in "ACGT" for base in fragment.sequence) for fragment in materialized
    )
    return {
        "records": len(materialized),
        "bases": total_bases,
        "gc_fraction": gc_bases / total_bases if total_bases else None,
        "ambiguous_fraction": ambiguous_bases / total_bases if total_bases else None,
        "records_by_label": dict(sorted(by_label.items())),
        "records_by_label_and_length": dict(sorted(by_label_and_length.items())),
        "records_by_length": dict(sorted(by_length.items(), key=lambda item: int(item[0]))),
        "source_genomes": len(by_genome),
        "records_by_genome": dict(sorted(by_genome.items())),
    }


_DEFAULT_CHECKSUM_EXCLUSIONS = ("checksums.sha256", "execution.json")
_VOLATILE_EXECUTION_RECORD = "execution.json"


def _resolved_bundle_root(root: Path) -> Path:
    """Return a real bundle directory without accepting a symlink as its root."""

    expanded = root.expanduser()
    if expanded.is_symlink():
        raise IntegrityError(f"Bundle root must not be a symbolic link: {expanded}")
    try:
        resolved = expanded.resolve(strict=True)
    except OSError as error:
        raise IntegrityError(f"Cannot access bundle root {expanded}: {error}") from error
    if not resolved.is_dir():
        raise IntegrityError(f"Bundle root is not a directory: {resolved}")
    return resolved


def _checksum_location(root: Path, requested: Path | None) -> tuple[Path, str]:
    """Resolve a checksum-manifest location confined to *root*."""

    checksum_path = requested.expanduser() if requested is not None else root / "checksums.sha256"
    if not checksum_path.is_absolute():
        checksum_path = Path.cwd() / checksum_path
    checksum_path = checksum_path.absolute()
    try:
        relative = checksum_path.relative_to(root).as_posix()
    except ValueError as error:
        raise IntegrityError(f"Checksum manifest must be inside bundle root {root}") from error
    _validate_relative_checksum_path(relative, path=checksum_path, line_number=None)
    return checksum_path, relative


def _validate_relative_checksum_path(
    relative: str,
    *,
    path: Path,
    line_number: int | None,
) -> str:
    """Validate one canonical, relative POSIX checksum path."""

    location = str(path) if line_number is None else f"{path}:{line_number}"
    if not relative:
        raise IntegrityError(f"{location}: checksum path must not be empty")
    if PurePosixPath(relative).is_absolute() or PureWindowsPath(relative).is_absolute():
        raise IntegrityError(f"{location}: absolute checksum paths are forbidden")
    if relative == ".." or relative.startswith("../"):
        raise IntegrityError(f"{location}: path escapes bundle root")
    posix = PurePosixPath(relative)
    if (
        "\\" in relative
        or any(character in relative for character in ("\x00", "\n", "\r"))
        or ".." in posix.parts
        or posix.as_posix() != relative
    ):
        raise IntegrityError(
            f"{location}: checksum path is not canonical POSIX syntax: {relative!r}"
        )
    return relative


def _bundle_file_inventory(
    root: Path,
    *,
    checksum_relative: str,
) -> dict[str, Path]:
    """Return the exact hashable inventory, rejecting every bundle symlink."""

    excluded = {checksum_relative, _VOLATILE_EXECUTION_RECORD}
    try:
        paths = sorted(root.rglob("*"), key=lambda item: item.as_posix())
    except OSError as error:
        raise IntegrityError(f"Cannot enumerate bundle {root}: {error}") from error

    inventory: dict[str, Path] = {}
    for path in paths:
        relative = path.relative_to(root).as_posix()
        _validate_relative_checksum_path(relative, path=path, line_number=None)
        if path.is_symlink():
            raise IntegrityError(f"Bundle contains a forbidden symbolic link: {relative!r}")
        if path.is_file() and relative not in excluded:
            inventory[relative] = path
    if not inventory:
        raise IntegrityError(
            "Checksum inventory is empty; a bundle must contain at least one "
            "regular file besides the checksum manifest and execution.json"
        )
    return inventory


def _validate_checksum_exclusions(exclude: Iterable[str], checksum_relative: str) -> None:
    """Retain the legacy writer argument without permitting integrity bypasses."""

    supplied = tuple(exclude)
    if len(supplied) != len(set(supplied)):
        raise ValueError("checksum exclusions must not contain duplicates")
    accepted = {
        frozenset(_DEFAULT_CHECKSUM_EXCLUSIONS),
        frozenset((checksum_relative, _VOLATILE_EXECUTION_RECORD)),
    }
    if frozenset(supplied) not in accepted:
        raise ValueError(
            "Only the checksum manifest and execution.json may be excluded from checksums"
        )


def write_checksums(
    root: Path,
    *,
    destination: Path | None = None,
    exclude: Iterable[str] = _DEFAULT_CHECKSUM_EXCLUSIONS,
) -> Path:
    """Write the exact, nonempty regular-file inventory in canonical order.

    The checksum manifest itself and the volatile ``execution.json`` record are
    the only permitted exclusions.  Symbolic links anywhere below *root* are
    rejected instead of silently omitted.
    """

    bundle_root = _resolved_bundle_root(root)
    output, checksum_relative = _checksum_location(bundle_root, destination)
    _validate_checksum_exclusions(exclude, checksum_relative)
    inventory = _bundle_file_inventory(bundle_root, checksum_relative=checksum_relative)
    rows = [f"{sha256_file(path)}  {relative}" for relative, path in inventory.items()]
    return write_text(output, "\n".join(rows))


def verify_checksums(root: Path, *, manifest: Path | None = None) -> dict[str, object]:
    """Verify a nonempty checksum manifest against the exact bundle inventory."""

    bundle_root = _resolved_bundle_root(root)
    checksum_path, checksum_relative = _checksum_location(bundle_root, manifest)
    inventory = _bundle_file_inventory(bundle_root, checksum_relative=checksum_relative)
    try:
        lines = checksum_path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise IntegrityError(f"Cannot read checksum manifest {checksum_path}: {error}") from error
    if not lines:
        raise IntegrityError(f"Checksum manifest {checksum_path} is empty")

    entries: dict[str, str] = {}
    listed_paths: list[str] = []
    for line_number, line in enumerate(lines, start=1):
        if not line:
            raise IntegrityError(f"{checksum_path}:{line_number}: blank checksum row")
        try:
            expected, relative = line.split("  ", 1)
        except ValueError as error:
            raise IntegrityError(
                f"{checksum_path}:{line_number}: malformed checksum row"
            ) from error
        if len(expected) != 64 or any(
            character not in "0123456789abcdef" for character in expected
        ):
            raise IntegrityError(f"{checksum_path}:{line_number}: invalid SHA-256 digest")
        relative = _validate_relative_checksum_path(
            relative,
            path=checksum_path,
            line_number=line_number,
        )
        if relative in entries:
            raise IntegrityError(f"{checksum_path}:{line_number}: duplicate path {relative!r}")
        entries[relative] = expected
        listed_paths.append(relative)

    if listed_paths != sorted(listed_paths):
        raise IntegrityError(f"{checksum_path}: checksum paths are not in canonical sorted order")

    expected_paths = set(inventory)
    recorded_paths = set(entries)
    missing = sorted(expected_paths - recorded_paths)
    extra = sorted(recorded_paths - expected_paths)
    if missing or extra:
        details: list[str] = []
        if missing:
            details.append("omitted file(s): " + ", ".join(missing[:10]))
        if extra:
            details.append("unexpected path(s): " + ", ".join(extra[:10]))
        raise IntegrityError("Checksum inventory mismatch; " + "; ".join(details))

    failures = [
        relative for relative, path in inventory.items() if sha256_file(path) != entries[relative]
    ]
    if failures:
        raise IntegrityError("Checksum verification failed for: " + ", ".join(failures[:10]))
    return {"status": "pass", "files_checked": len(inventory)}


def summarize_split_truth(rows: Iterable[Mapping[str, object]]) -> dict[str, object]:
    """Summarize already-serialized truth without re-reading FASTA."""

    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        counts[str(row["partition"])][str(row["label"])] += 1
    return {
        partition: dict(sorted(label_counts.items()))
        for partition, label_counts in sorted(counts.items())
    }
