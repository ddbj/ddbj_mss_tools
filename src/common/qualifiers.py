"""
Source-feature qualifier definitions shared by DDBJ MSS tools.

Loads source_qualifiers.json (a copy of the GGS web frontend's list) and
exposes lookups by qualifier name. Each definition is a dict with at least
``name``, ``boolean`` (flag-type qualifier without a value) and ``multiple``
(may appear more than once).
"""

from __future__ import annotations

import json
from pathlib import Path

_DEFINITIONS_PATH = Path(__file__).parent / "data" / "source_qualifiers.json"
_by_name: dict[str, dict] | None = None


def _load() -> dict[str, dict]:
    global _by_name
    if _by_name is None:
        items = json.loads(_DEFINITIONS_PATH.read_text(encoding="utf-8"))
        _by_name = {item["name"]: item for item in items}
    return _by_name


def known_qualifier_names() -> frozenset[str]:
    """Return the names of all known source qualifiers."""
    return frozenset(_load())


def qualifier_def(name: str) -> dict | None:
    """Return the definition of qualifier *name*, or None if it is unknown."""
    return _load().get(name)
