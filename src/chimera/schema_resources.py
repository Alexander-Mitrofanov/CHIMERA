"""Versioned JSON Schema discovery and offline instance validation.

Schemas are loaded from package data in an installed wheel.  Editable source
checkouts fall back to the repository's top-level ``schemas`` directory, so the
same public API works in both development and publication environments.
"""

from __future__ import annotations

import json
from enum import StrEnum
from importlib import resources
from importlib.resources.abc import Traversable
from pathlib import Path
from typing import Final, TypeAlias, TypeGuard, cast

from jsonschema import Draft202012Validator, FormatChecker
from jsonschema.exceptions import SchemaError
from referencing import Registry, Resource
from referencing.jsonschema import DRAFT202012, SchemaRegistry

from .errors import ConfigurationError, IntegrityError

JsonScalar: TypeAlias = str | int | float | bool | None
JsonValue: TypeAlias = JsonScalar | list["JsonValue"] | dict[str, "JsonValue"]
JsonObject: TypeAlias = dict[str, JsonValue]


class SchemaName(StrEnum):
    """Stable public names for CHIMERA's versioned JSON Schemas."""

    BUNDLE = "bundle"
    SPLIT = "split"
    METADATA_ROW = "metadata-row"
    TRUTH_ROW = "truth-row"
    REFERENCE_ROW = "reference-row"
    SEQUENCE_ROW = "sequence-row"
    ASSIGNMENT_ROW = "assignment-row"
    EXCLUSION_ROW = "exclusion-row"
    RESOLVED_CONFIG = "resolved-config"


_SCHEMA_FILES: Final[dict[SchemaName, str]] = {
    SchemaName.BUNDLE: "benchmark-bundle.schema.json",
    SchemaName.SPLIT: "split-manifest.schema.json",
    SchemaName.METADATA_ROW: "reference-metadata-row.schema.json",
    SchemaName.TRUTH_ROW: "truth-row.schema.json",
    SchemaName.REFERENCE_ROW: "reference-row.schema.json",
    SchemaName.SEQUENCE_ROW: "sequence-row.schema.json",
    SchemaName.ASSIGNMENT_ROW: "assignment-row.schema.json",
    SchemaName.EXCLUSION_ROW: "exclusion-row.schema.json",
    SchemaName.RESOLVED_CONFIG: "resolved-config.schema.json",
}

JSON_SCHEMA_NAMES: Final[tuple[str, ...]] = tuple(name.value for name in SchemaName)
"""Canonical names accepted by ``chimera schema`` and :func:`load_schema`."""


class SchemaValidationError(IntegrityError):
    """A JSON-compatible instance violates a published CHIMERA schema."""


def _coerce_schema_name(name: SchemaName | str) -> SchemaName:
    if isinstance(name, SchemaName):
        return name
    try:
        return SchemaName(name)
    except (TypeError, ValueError) as error:
        choices = ", ".join(JSON_SCHEMA_NAMES)
        raise ConfigurationError(f"Unknown JSON Schema {name!r}; choose from: {choices}") from error


def _schema_candidates(filename: str) -> tuple[Traversable, ...]:
    packaged = resources.files("chimera").joinpath("schemas", filename)
    source_checkout = Path(__file__).resolve().parents[2] / "schemas" / filename
    return packaged, source_checkout


def _read_schema_text(name: SchemaName) -> str:
    filename = _SCHEMA_FILES[name]
    failures: list[str] = []
    for candidate in _schema_candidates(filename):
        try:
            if candidate.is_file():
                return candidate.read_text(encoding="utf-8")
        except (OSError, UnicodeError) as error:
            failures.append(f"{candidate}: {error}")
    detail = f" ({'; '.join(failures)})" if failures else ""
    raise ConfigurationError(f"Packaged JSON Schema resource {filename!r} is unavailable{detail}")


def _is_json_value(value: object) -> TypeGuard[JsonValue]:
    if value is None or isinstance(value, (str, int, float, bool)):
        return True
    if isinstance(value, list):
        return all(_is_json_value(item) for item in value)
    if isinstance(value, dict):
        return all(isinstance(key, str) and _is_json_value(item) for key, item in value.items())
    return False


def load_schema(name: SchemaName | str) -> JsonObject:
    """Load one schema as a fresh JSON object.

    Returning a fresh object prevents mutations by one caller from corrupting
    subsequent validation or CLI output.
    """

    schema_name = _coerce_schema_name(name)
    try:
        decoded: object = json.loads(_read_schema_text(schema_name))
    except json.JSONDecodeError as error:
        raise ConfigurationError(
            f"Packaged JSON Schema {_SCHEMA_FILES[schema_name]!r} is malformed: {error}"
        ) from error
    if not isinstance(decoded, dict) or not _is_json_value(decoded):
        raise ConfigurationError(
            f"Packaged JSON Schema {_SCHEMA_FILES[schema_name]!r} is not a JSON object"
        )
    return cast(JsonObject, decoded)


def schema_filename(name: SchemaName | str) -> str:
    """Return the canonical package filename for a public schema name."""

    return _SCHEMA_FILES[_coerce_schema_name(name)]


def _offline_registry() -> SchemaRegistry:
    registry: SchemaRegistry = Registry()
    for name in SchemaName:
        document = load_schema(name)
        identifier = document.get("$id")
        if not isinstance(identifier, str) or not identifier:
            raise ConfigurationError(f"JSON Schema {name.value!r} has no nonempty $id")
        registry = registry.with_resource(
            identifier,
            Resource.from_contents(document, default_specification=DRAFT202012),
        )
    return registry


def validate_instance(instance: object, schema: SchemaName | str) -> None:
    """Validate *instance* against one packaged schema without network access.

    Raises:
        ConfigurationError: If a packaged schema is missing or invalid.
        SchemaValidationError: If *instance* violates the selected contract.
    """

    schema_name = _coerce_schema_name(schema)
    document = load_schema(schema_name)
    try:
        Draft202012Validator.check_schema(document)
    except SchemaError as error:
        raise ConfigurationError(
            f"Packaged JSON Schema {schema_name.value!r} is invalid: {error.message}"
        ) from error
    validator = Draft202012Validator(
        document,
        registry=_offline_registry(),
        format_checker=FormatChecker(),
    )
    errors = sorted(
        validator.iter_errors(instance),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if errors:
        validation_error = errors[0]
        location = "/".join(str(part) for part in validation_error.absolute_path) or "$"
        raise SchemaValidationError(
            f"{schema_name.value} schema violation at {location}: {validation_error.message}"
        )


def validate_packaged_schemas() -> None:
    """Check every packaged schema and all in-package references."""

    registry = _offline_registry()
    for name in SchemaName:
        document = load_schema(name)
        try:
            Draft202012Validator.check_schema(document)
            Draft202012Validator(document, registry=registry)
        except SchemaError as error:
            raise ConfigurationError(
                f"Packaged JSON Schema {name.value!r} is invalid: {error.message}"
            ) from error


__all__ = [
    "JSON_SCHEMA_NAMES",
    "JsonObject",
    "JsonValue",
    "SchemaName",
    "SchemaValidationError",
    "load_schema",
    "schema_filename",
    "validate_instance",
    "validate_packaged_schemas",
]
