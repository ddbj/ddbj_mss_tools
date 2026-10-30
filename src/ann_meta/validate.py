"""Input checks for `ann_meta apply`: the common JSON and the entry table.

``validate`` never raises on bad input; it returns every problem as an
:class:`Issue` whose ``message`` is a short machine-friendly code:

* ``forbidden_char`` — a control character (Unicode category ``Cc``: tab, CR,
  LF, ``\\x7f-\\x9f`` ...) or anything ``str.splitlines`` treats as a line
  break (``\\x1c-\\x1e``, ``\\u2028``, ``\\u2029`` ...) in any value or
  qualifier name. ``AnnFile.read`` splits lines with ``splitlines``, so such a
  value could move rows into another feature or entry.
* ``unknown_qualifier`` — a qualifier name in ``SOURCE``, an entry's overrides
  or an entry's ``extra`` that is neither in ``source_qualifiers.json`` nor in
  *allowed_unknown* (``MetaDoc.unknown_qualifiers`` on a round trip).
* ``reserved_qualifier`` — ``ff_definition`` / ``submitter_seqid`` / a role name
  given in ``SOURCE`` or overrides (apply generates these itself).
* ``bad_type`` / ``bad_status`` / ``bad_topology`` — value outside the domain.
* role rules — ``unplaced_with_name``, ``unplaced_not_partial``,
  ``unplaced_circular``, ``partial_circular``, ``plasmid_needs_name``,
  ``organelle_needs_name``, ``chromosome_needs_name`` / ``segment_needs_name``
  (only when there are two or more of them).
* ``role_not_allowed_for_category`` — a role other than ``unplaced`` for a
  category whose entries are all unplaced (WGS family, ENV, TSA, GENE).
* ``empty_seq_id``, ``duplicate_seq_id``, ``seq_id_not_in_fasta``, ``no_entries``,
  ``seq_id_missing_from_table`` (a FASTA id with no row in the table).
* ``unknown_category`` — not ``NONE``, ``GENE`` or a name in
  ``submission_categories.json``.
"""
from __future__ import annotations

import unicodedata
from dataclasses import dataclass

from ann_meta.categories import PASS_THROUGH, ROLE_CATEGORIES
from ann_meta.model import Entry
from common.qualifiers import known_qualifier_names
from common.submission_category import get_category_rules

ROLES = ("chromosome", "plasmid", "organelle", "segment")
TYPES = ROLES + ("unplaced",)
STATUSES = ("complete", "partial")
TOPOLOGIES = ("circular", "linear", "")
RESERVED = frozenset({"ff_definition", "submitter_seqid", *ROLES})


@dataclass
class Issue:
    field: str
    entry: str | None
    message: str


def _bad(value) -> bool:
    if not isinstance(value, str):
        return False
    if len(f"x{value}x".splitlines()) > 1:
        return True
    return any(unicodedata.category(c) == "Cc" for c in value)


def _walk(value, path: str, issues: list[Issue]) -> None:
    """Report forbidden characters in every string (and every dict key) under *value*."""
    if isinstance(value, dict):
        for k, v in value.items():
            sub = f"{path}.{k}" if path else str(k)
            if _bad(k):
                issues.append(Issue(sub, None, "forbidden_char"))
            _walk(v, sub, issues)
    elif isinstance(value, (list, tuple)):
        for i, v in enumerate(value):
            _walk(v, f"{path}.{i}", issues)
    elif _bad(value):
        issues.append(Issue(path, None, "forbidden_char"))


def _category_ok(category: str) -> bool:
    if category in PASS_THROUGH:
        return True
    try:
        get_category_rules(category)
    except ValueError:
        return False
    return True


def validate(common: dict, entries: list[Entry], *, category: str, fasta_ids: list[str] | None,
             allowed_unknown: set[str]) -> list[Issue]:
    issues: list[Issue] = []
    if not _category_ok(category):
        issues.append(Issue("category", None, "unknown_category"))

    _walk(common, "", issues)
    known = known_qualifier_names()
    allowed = set(known) | set(allowed_unknown)

    source = common.get("SOURCE") or {}
    if isinstance(source, dict):
        for q in source:
            if q in RESERVED:
                issues.append(Issue(f"SOURCE.{q}", None, "reserved_qualifier"))
            elif q not in allowed:
                issues.append(Issue(f"SOURCE.{q}", None, "unknown_qualifier"))

    if not entries:
        issues.append(Issue("entries", None, "no_entries"))
        return issues

    counts = {r: sum(1 for e in entries if e.type == r) for r in ROLES}
    seen: set[str] = set()
    fasta = set(fasta_ids) if fasta_ids is not None else None
    for e in entries:
        eid = e.seq_id or None
        cells = [("seq_id", e.seq_id), ("type", e.type), ("seq_name", e.seq_name), ("status", e.status),
                 ("topology", e.topology), ("ff_definition", e.ff_definition)]
        for q, v in [*e.overrides.items(), *e.extra]:
            if _bad(q):
                issues.append(Issue(q, eid, "forbidden_char"))
            cells.append((q, v))
        for name, v in cells:
            if _bad(v):
                issues.append(Issue(name, eid, "forbidden_char"))
        for q in e.overrides:
            if q in RESERVED:
                issues.append(Issue(q, eid, "reserved_qualifier"))
            elif q not in allowed:
                issues.append(Issue(q, eid, "unknown_qualifier"))
        for q in dict.fromkeys(q for q, _ in e.extra):
            if q not in allowed:
                issues.append(Issue(q, eid, "unknown_qualifier"))

        if not e.seq_id:
            issues.append(Issue("seq_id", None, "empty_seq_id"))
        elif e.seq_id in seen:
            issues.append(Issue("seq_id", eid, "duplicate_seq_id"))
        else:
            seen.add(e.seq_id)
            if fasta is not None and e.seq_id not in fasta:
                issues.append(Issue("seq_id", eid, "seq_id_not_in_fasta"))

        if e.type not in TYPES:
            issues.append(Issue("type", eid, "bad_type"))
        if e.status not in STATUSES:
            issues.append(Issue("status", eid, "bad_status"))
        if e.topology not in TOPOLOGIES:
            issues.append(Issue("topology", eid, "bad_topology"))

        if e.type == "unplaced":
            if e.seq_name:
                issues.append(Issue("seq_name", eid, "unplaced_with_name"))
            if e.status != "partial":
                issues.append(Issue("status", eid, "unplaced_not_partial"))
            if e.topology == "circular":
                issues.append(Issue("topology", eid, "unplaced_circular"))
        else:
            if e.status == "partial" and e.topology == "circular":
                issues.append(Issue("topology", eid, "partial_circular"))
            if e.type in ("plasmid", "organelle") and not e.seq_name:
                issues.append(Issue("seq_name", eid, f"{e.type}_needs_name"))
            if e.type in ("chromosome", "segment") and counts[e.type] >= 2 and not e.seq_name:
                issues.append(Issue("seq_name", eid, f"{e.type}_needs_name"))
            if e.type in ROLES and category not in ROLE_CATEGORIES:
                issues.append(Issue("type", eid, "role_not_allowed_for_category"))
    if fasta_ids is not None:
        for fid in fasta_ids:
            if fid not in seen:
                issues.append(Issue("entries", fid, "seq_id_missing_from_table"))
    return issues
