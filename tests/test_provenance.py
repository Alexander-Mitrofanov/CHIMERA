from __future__ import annotations

import re

from chimera.provenance import software_content_sha256


def test_software_content_receipt_is_stable_sha256() -> None:
    first = software_content_sha256()
    second = software_content_sha256()

    assert first == second
    assert re.fullmatch(r"[0-9a-f]{64}", first)
