"""Validated configuration for reproducible benchmark generation.

The command line and TOML configuration file both resolve into the same immutable
``BenchmarkConfig``.  Keeping scientific defaults here prevents the CLI and library
API from silently diverging.
"""

from __future__ import annotations

import math
import tomllib
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, field
from datetime import date
from enum import StrEnum
from pathlib import Path
from typing import Any, Literal, Self

from .errors import ConfigurationError

CONFIG_SCHEMA_VERSION = 1


class SplitKind(StrEnum):
    """Supported evaluation protocols, named after the study design."""

    RANDOM = "random"
    GENOME = "genome"
    SIMILARITY = "similarity"
    TEMPORAL = "temporal"
    TAXONOMY = "taxonomy"

    @property
    def directory_name(self) -> str:
        """Return the stable, human-sortable output directory name."""

        prefixes = {
            SplitKind.RANDOM: "2a_random_fragment",
            SplitKind.GENOME: "2b_genome_holdout",
            SplitKind.SIMILARITY: "2c_similarity_filtered",
            SplitKind.TEMPORAL: "2d_temporal_holdout",
            SplitKind.TAXONOMY: "2e_taxonomic_holdout",
        }
        return prefixes[self]

    @classmethod
    def parse(cls, value: str) -> Self:
        """Parse friendly names and protocol identifiers."""

        normalized = value.strip().lower().replace("_", "-")
        aliases = {
            "2a": cls.RANDOM,
            "random": cls.RANDOM,
            "random-fragment": cls.RANDOM,
            "2b": cls.GENOME,
            "genome": cls.GENOME,
            "genome-level": cls.GENOME,
            "2c": cls.SIMILARITY,
            "similarity": cls.SIMILARITY,
            "similarity-filtered": cls.SIMILARITY,
            "2d": cls.TEMPORAL,
            "temporal": cls.TEMPORAL,
            "2e": cls.TAXONOMY,
            "taxonomy": cls.TAXONOMY,
            "taxonomic": cls.TAXONOMY,
            "taxonomic-holdout": cls.TAXONOMY,
        }
        try:
            return aliases[normalized]
        except KeyError as exc:
            choices = ", ".join(kind.value for kind in cls)
            raise ConfigurationError(
                f"Unknown split {value!r}; choose from: {choices}, all"
            ) from exc


class DuplicatePolicy(StrEnum):
    """Policy for duplicate genomes within one class."""

    ERROR = "error"
    DROP = "drop"


class MissingMetadataPolicy(StrEnum):
    """Policy for records missing dates or taxonomy needed by a split."""

    ERROR = "error"
    EXCLUDE = "exclude"


@dataclass(frozen=True, slots=True)
class SimilarityBands:
    """Novelty bands expressed as estimated nucleotide similarity in [0, 1]."""

    high: float = 0.90
    moderate: float = 0.70
    low: float = 0.30

    def __post_init__(self) -> None:
        values = (self.high, self.moderate, self.low)
        if any(not math.isfinite(value) for value in values):
            raise ConfigurationError("Similarity band thresholds must be finite")
        if not 1.0 >= self.high > self.moderate > self.low >= 0.0:
            raise ConfigurationError(
                "Similarity bands must satisfy 1 >= high > moderate > low >= 0"
            )

    def classify(self, similarity: float | None) -> str:
        """Map a maximum train similarity to one mutually exclusive similarity bin."""

        if similarity is None:
            return "no_detectable_match"
        if similarity < self.low:
            return "distant_detectable"
        if similarity < self.moderate:
            return "low_similarity"
        if similarity < self.high:
            return "moderate_similarity"
        return "high_similarity"


def _parse_date(value: object, field_name: str) -> date | None:
    if value is None or value == "":
        return None
    if isinstance(value, date):
        return value
    try:
        return date.fromisoformat(str(value))
    except ValueError as exc:
        raise ConfigurationError(f"{field_name} must be an ISO date (YYYY-MM-DD)") from exc


def _as_paths(value: object, *, base_dir: Path, field_name: str) -> tuple[Path, ...]:
    if value is None:
        return ()
    if isinstance(value, (str, Path)):
        values = [value]
    elif isinstance(value, Iterable) and not isinstance(value, Mapping):
        values = list(value)
    else:
        raise ConfigurationError(f"{field_name} must be a path or a list of paths")
    paths: list[Path] = []
    for item in values:
        path = Path(str(item)).expanduser()
        if not path.is_absolute():
            path = base_dir / path
        paths.append(path.resolve())
    if len(set(paths)) != len(paths):
        raise ConfigurationError(f"{field_name} contains the same path more than once")
    return tuple(paths)


def _as_output_path(value: object, *, base_dir: Path) -> Path:
    """Parse exactly one non-empty output directory path."""

    if isinstance(value, (str, Path)):
        items = (value,)
    elif isinstance(value, Iterable) and not isinstance(value, (bytes, Mapping)):
        items = tuple(value)
    else:
        raise ConfigurationError("output_dir must be exactly one non-empty path")
    if len(items) != 1:
        raise ConfigurationError("output_dir must be exactly one non-empty path")
    item = items[0]
    if (
        not isinstance(item, (str, Path))
        or (isinstance(item, str) and not item.strip())
        or Path(item) == Path()
    ):
        raise ConfigurationError("output_dir must be exactly one non-empty path")
    return _as_paths(item, base_dir=base_dir, field_name="output_dir")[0]


def _parse_fragment_lengths(value: object) -> tuple[int, ...]:
    """Parse one fragment length or an iterable of fragment lengths."""

    if isinstance(value, int) and not isinstance(value, bool):
        values = (value,)
    elif isinstance(value, Iterable) and not isinstance(value, (str, bytes, Mapping)):
        values = tuple(value)
    else:
        raise ConfigurationError("fragment_lengths must be an integer or a list of integers")
    try:
        return tuple(int(item) for item in values)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(
            "fragment_lengths must be an integer or a list of integers"
        ) from exc


def _parse_holdout_taxa(value: object) -> tuple[str, ...]:
    """Parse one taxon name or an iterable of taxon names."""

    if isinstance(value, str):
        values = (value,)
    elif isinstance(value, Iterable) and not isinstance(value, (bytes, Mapping)):
        values = tuple(value)
    else:
        raise ConfigurationError("holdout_taxa must be a name or a list of names")
    return tuple(str(item) for item in values)


def _parse_similarity_bands(value: object) -> SimilarityBands:
    """Parse the optional similarity-band table without leaking ``TypeError``."""

    if value is None:
        return SimilarityBands()
    if not isinstance(value, Mapping):
        raise ConfigurationError(
            "similarity_bands must be a table containing numeric high, moderate, and low values"
        )
    keys = set(value)
    allowed = {"high", "moderate", "low"}
    if any(not isinstance(key, str) for key in keys) or not keys <= allowed:
        raise ConfigurationError(
            "similarity_bands may contain only high, moderate, and low thresholds"
        )
    try:
        return SimilarityBands(**dict(value))
    except TypeError as exc:
        raise ConfigurationError(
            "similarity_bands must contain numeric high, moderate, and low values"
        ) from exc


def parse_splits(value: object) -> tuple[SplitKind, ...]:
    """Parse a split list while preserving canonical protocol order."""

    if value is None or value == "all" or value == ["all"]:
        return tuple(SplitKind)
    if isinstance(value, str):
        raw = [value]
    elif isinstance(value, Iterable) and not isinstance(value, Mapping):
        raw = list(value)
    else:
        raise ConfigurationError("splits must be a name or a list of names")
    parsed: set[SplitKind] = set()
    for item in raw:
        for token in str(item).split(","):
            if token.strip():
                parsed.add(SplitKind.parse(token))
    if not parsed:
        raise ConfigurationError("At least one split must be selected")
    return tuple(kind for kind in SplitKind if kind in parsed)


@dataclass(frozen=True, slots=True)
class BenchmarkConfig:
    """Complete, validated parameters for one benchmark bundle."""

    virus_paths: tuple[Path, ...]
    host_paths: tuple[Path, ...]
    output_dir: Path
    metadata_path: Path | None = None
    splits: tuple[SplitKind, ...] = field(default_factory=lambda: tuple(SplitKind))
    seed: int = 42
    test_fraction: float = 0.20
    fragment_lengths: tuple[int, ...] = (500,)
    fragments_per_genome: int = 100
    strand_mode: Literal["both", "forward"] = "both"
    max_ambiguous_fraction: float = 0.05
    duplicate_policy: DuplicatePolicy = DuplicatePolicy.ERROR
    missing_metadata: MissingMetadataPolicy = MissingMetadataPolicy.ERROR
    temporal_cutoff: date | None = None
    taxonomy_rank: str = "family"
    holdout_taxa: tuple[str, ...] = ()
    auto_holdout_count: int = 1
    similarity_k: int = 21
    sketch_size: int = 2_000
    max_train_similarity: float = 0.95
    min_similarity_coverage: float = 0.85
    similarity_bands: SimilarityBands = field(default_factory=SimilarityBands)
    similarity_table: Path | None = None
    overwrite: bool = False

    def __post_init__(self) -> None:
        if not self.virus_paths:
            raise ConfigurationError("At least one --virus FASTA file or directory is required")
        if not self.host_paths:
            raise ConfigurationError("At least one --host FASTA file or directory is required")
        if not self.splits:
            raise ConfigurationError("At least one evaluation split is required")
        if not 0.0 < self.test_fraction < 1.0:
            raise ConfigurationError("test_fraction must be greater than 0 and less than 1")
        if not self.fragment_lengths or any(length < 1 for length in self.fragment_lengths):
            raise ConfigurationError("fragment_lengths must contain positive integers")
        if len(set(self.fragment_lengths)) != len(self.fragment_lengths):
            raise ConfigurationError("fragment_lengths must not contain duplicates")
        minimum_fragments = 2 * len(self.fragment_lengths)
        if self.fragments_per_genome < minimum_fragments:
            raise ConfigurationError(
                f"fragments_per_genome must be at least {minimum_fragments} so every "
                "genome contributes train and test records at every requested length"
            )
        if self.strand_mode not in {"both", "forward"}:
            raise ConfigurationError("strand_mode must be 'both' or 'forward'")
        if not 0.0 <= self.max_ambiguous_fraction <= 1.0:
            raise ConfigurationError("max_ambiguous_fraction must be between 0 and 1")
        if not self.taxonomy_rank or any(char.isspace() for char in self.taxonomy_rank):
            raise ConfigurationError("taxonomy_rank must be one non-empty field name")
        if self.auto_holdout_count < 1:
            raise ConfigurationError("auto_holdout_count must be at least 1")
        if not 5 <= self.similarity_k <= 63 or self.similarity_k % 2 == 0:
            raise ConfigurationError("similarity_k must be an odd integer between 5 and 63")
        if self.sketch_size < 100:
            raise ConfigurationError("sketch_size must be at least 100")
        if not 0.0 <= self.max_train_similarity <= 1.0:
            raise ConfigurationError("max_train_similarity must be between 0 and 1")
        if not 0.0 <= self.min_similarity_coverage <= 1.0:
            raise ConfigurationError("min_similarity_coverage must be between 0 and 1")
        if self.max_train_similarity < self.similarity_bands.high:
            raise ConfigurationError(
                "max_train_similarity must be at least the high similarity-band threshold"
            )
        normalized_rank = self.taxonomy_rank.strip().lower()
        object.__setattr__(self, "taxonomy_rank", normalized_rank)
        object.__setattr__(self, "fragment_lengths", tuple(sorted(self.fragment_lengths)))
        object.__setattr__(
            self,
            "holdout_taxa",
            tuple(dict.fromkeys(taxon.strip() for taxon in self.holdout_taxa if taxon.strip())),
        )

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any], *, base_dir: Path | None = None) -> Self:
        """Construct a config from TOML/JSON-like values.

        Hyphenated CLI-style keys are accepted and converted to underscores.
        Relative paths are resolved from the configuration file's directory.
        """

        raw = {str(key).replace("-", "_"): value for key, value in values.items()}
        known = {field_.name for field_ in cls.__dataclass_fields__.values()}
        unknown = sorted(set(raw) - known - {"schema_version"})
        if unknown:
            raise ConfigurationError(f"Unknown configuration field(s): {', '.join(unknown)}")
        schema_version = int(raw.pop("schema_version", CONFIG_SCHEMA_VERSION))
        if schema_version != CONFIG_SCHEMA_VERSION:
            raise ConfigurationError(
                f"Unsupported config schema_version {schema_version}; expected {CONFIG_SCHEMA_VERSION}"
            )
        path_base = (base_dir or Path.cwd()).resolve()
        virus_paths = _as_paths(
            raw.pop("virus_paths", ()), base_dir=path_base, field_name="virus_paths"
        )
        host_paths = _as_paths(
            raw.pop("host_paths", ()), base_dir=path_base, field_name="host_paths"
        )
        output_raw = raw.pop("output_dir", "benchmark-output")
        output_dir = _as_output_path(output_raw, base_dir=path_base)
        metadata_values = _as_paths(
            raw.pop("metadata_path", None), base_dir=path_base, field_name="metadata_path"
        )
        similarity_values = _as_paths(
            raw.pop("similarity_table", None), base_dir=path_base, field_name="similarity_table"
        )
        bands_raw = raw.pop("similarity_bands", None)
        bands = _parse_similarity_bands(bands_raw)
        duplicate = DuplicatePolicy(str(raw.pop("duplicate_policy", DuplicatePolicy.ERROR.value)))
        missing = MissingMetadataPolicy(
            str(raw.pop("missing_metadata", MissingMetadataPolicy.ERROR.value))
        )
        fragment_lengths = _parse_fragment_lengths(raw.pop("fragment_lengths", (500,)))
        holdout_taxa = _parse_holdout_taxa(raw.pop("holdout_taxa", ()))
        return cls(
            virus_paths=virus_paths,
            host_paths=host_paths,
            output_dir=output_dir,
            metadata_path=metadata_values[0] if metadata_values else None,
            splits=parse_splits(raw.pop("splits", None)),
            temporal_cutoff=_parse_date(raw.pop("temporal_cutoff", None), "temporal_cutoff"),
            fragment_lengths=fragment_lengths,
            holdout_taxa=holdout_taxa,
            duplicate_policy=duplicate,
            missing_metadata=missing,
            similarity_bands=bands,
            similarity_table=similarity_values[0] if similarity_values else None,
            **raw,
        )

    @classmethod
    def from_toml(cls, path: Path) -> Self:
        """Load ``[benchmark]`` from a TOML file with relative path support."""

        try:
            with path.open("rb") as handle:
                document = tomllib.load(handle)
        except (OSError, tomllib.TOMLDecodeError) as exc:
            raise ConfigurationError(f"Cannot read TOML configuration {path}: {exc}") from exc
        section = document.get("benchmark", document)
        if not isinstance(section, Mapping):
            raise ConfigurationError("TOML [benchmark] must be a table")
        return cls.from_mapping(section, base_dir=path.parent)

    def as_manifest_dict(self) -> dict[str, Any]:
        """Return a stable JSON-compatible representation for provenance."""

        result = asdict(self)
        result["schema_version"] = CONFIG_SCHEMA_VERSION
        for key in ("virus_paths", "host_paths"):
            result[key] = [str(path) for path in result[key]]
        for key in ("output_dir", "metadata_path", "similarity_table"):
            if result[key] is not None:
                result[key] = str(result[key])
        result["splits"] = [kind.value for kind in self.splits]
        result["temporal_cutoff"] = (
            self.temporal_cutoff.isoformat() if self.temporal_cutoff else None
        )
        result["duplicate_policy"] = self.duplicate_policy.value
        result["missing_metadata"] = self.missing_metadata.value
        return result
