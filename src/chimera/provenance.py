"""Content-addressed provenance for the installed CHIMERA implementation."""

from __future__ import annotations

import hashlib
import json
from importlib import resources

from .schema_resources import JSON_SCHEMA_NAMES, load_schema, schema_filename

_DOMAIN = b"CHIMERA-SOFTWARE-CONTENT-SHA256-v1\0"


def software_content_sha256() -> str:
    """Hash executable package sources and canonical schemas under logical paths.

    Unlike a Git revision, this receipt remains available in wheels and records
    the implementation bytes that actually generated or validated a bundle.
    """

    package_root = resources.files("chimera")
    payloads = {
        child.name: child.read_bytes()
        for child in package_root.iterdir()
        if child.is_file() and (child.name.endswith(".py") or child.name == "py.typed")
    }
    for schema_name in JSON_SCHEMA_NAMES:
        payloads[f"schemas/{schema_filename(schema_name)}"] = json.dumps(
            load_schema(schema_name),
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")

    digest = hashlib.sha256(_DOMAIN)
    for logical_path, payload in sorted(payloads.items()):
        path_bytes = logical_path.encode("utf-8")
        digest.update(len(path_bytes).to_bytes(4, "big"))
        digest.update(path_bytes)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()


__all__ = ["software_content_sha256"]
