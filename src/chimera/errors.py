"""Typed, actionable errors raised by CHIMERA."""

from __future__ import annotations


class ChimeraError(Exception):
    """Base class for expected user-facing CHIMERA failures."""


class ConfigurationError(ChimeraError):
    """A configuration is missing information or is scientifically invalid."""


class InputError(ChimeraError):
    """An input file or record violates the documented schema."""


class IntegrityError(ChimeraError):
    """A leakage, checksum, or output-integrity invariant failed."""


class ExternalToolError(ChimeraError):
    """An optional external bioinformatics tool could not be used safely."""
