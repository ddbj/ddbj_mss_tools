"""Parse MSS location strings and convert to BioPython SeqFeature locations."""

from __future__ import annotations

import re

from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    ExactPosition,
    SimpleLocation,
)


def expand_location(loc_str: str, seq_len: int, entry_id: str = "") -> str:
    """Replace meta-notations in location strings.

    E  → seq_len (as end position)
    @@[entry]@@ → entry_id
    """
    loc_str = re.sub(r"(?<!\w)E(?!\w)", str(seq_len), loc_str)
    if entry_id:
        loc_str = loc_str.replace("@@[entry]@@", entry_id)
    return loc_str


def parse_mss_location(loc_str: str) -> SimpleLocation | CompoundLocation:
    """Parse an MSS location string (1-based, inclusive) to a BioPython location.

    Handles: simple, complement, join, complement(join(...)), order.
    Partial positions: <start, >end.
    """
    return _parse_loc(loc_str.strip())


# ── Internal helpers ──────────────────────────────────────────────────────────

def _make_start(s: str) -> int | BeforePosition:
    s = s.strip()
    if s.startswith("<"):
        return BeforePosition(int(s[1:]) - 1)
    return ExactPosition(int(s) - 1)


def _make_end(s: str) -> int | AfterPosition:
    s = s.strip()
    if s.startswith(">"):
        return AfterPosition(int(s[1:]))
    return ExactPosition(int(s))


def _parse_simple_range(loc_str: str) -> tuple:
    m = re.match(r"([<>]?\d+)\.\.([<>]?\d+)", loc_str.strip())
    if not m:
        raise ValueError(f"Cannot parse location range: {loc_str!r}")
    return _make_start(m.group(1)), _make_end(m.group(2))


def _split_parts(inner: str) -> list[str]:
    """Split on commas that are not inside parentheses."""
    parts, depth, cur = [], 0, ""
    for ch in inner:
        if ch == "(":
            depth += 1
            cur += ch
        elif ch == ")":
            depth -= 1
            cur += ch
        elif ch == "," and depth == 0:
            parts.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        parts.append(cur.strip())
    return parts


def _parse_loc(loc_str: str) -> SimpleLocation | CompoundLocation:
    loc_str = loc_str.strip()

    if loc_str.startswith("complement(") and loc_str.endswith(")"):
        inner = loc_str[11:-1]
        if inner.startswith("join(") and inner.endswith(")"):
            sub_parts = [_parse_loc(p) for p in _split_parts(inner[5:-1])]
            locs = [SimpleLocation(p.start, p.end, -1) for p in sub_parts]
            return CompoundLocation(locs)
        else:
            start, end = _parse_simple_range(inner)
            return SimpleLocation(start, end, -1)

    if loc_str.startswith("join(") and loc_str.endswith(")"):
        sub_parts = [_parse_loc(p) for p in _split_parts(loc_str[5:-1])]
        return CompoundLocation(sub_parts)

    if loc_str.startswith("order(") and loc_str.endswith(")"):
        sub_parts = [_parse_loc(p) for p in _split_parts(loc_str[6:-1])]
        return CompoundLocation(sub_parts, operator="order")

    if ".." in loc_str:
        start, end = _parse_simple_range(loc_str)
        return SimpleLocation(start, end, +1)

    # Single position
    val = int(loc_str.strip("<>"))
    return SimpleLocation(val - 1, val, +1)
