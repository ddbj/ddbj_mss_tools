"""Input checks for `ann_meta apply`: the common JSON and the entry table.

``validate`` never raises on bad input; it returns every problem as an
:class:`Issue` whose ``message`` is a short machine-friendly code:

* ``forbidden_char`` — a tab, newline or other control character (``\\x00-\\x1f``,
  ``\\x7f``) in any value or qualifier name (they would break the 5-column .ann).
* ``unknown_qualifier`` — a qualifier name in ``SOURCE`` or an entry's overrides
  that is neither in ``source_qualifiers.json`` nor in *allowed_unknown*
  (names in an entry's ``extra`` are pass-through and always allowed).
* ``reserved_qualifier`` — ``ff_definition`` / ``submitter_seqid`` / a role name
  given in ``SOURCE`` or overrides (apply generates these itself).
* ``bad_type`` / ``bad_status`` / ``bad_topology`` — value outside the domain.
* role rules — ``unplaced_with_name``, ``unplaced_not_partial``,
  ``unplaced_circular``, ``partial_circular``, ``plasmid_needs_name``,
  ``organelle_needs_name``, ``chromosome_needs_name`` / ``segment_needs_name``
  (only when there are two or more of them).
* ``role_not_allowed_for_category`` — a role other than ``unplaced`` for a
  category whose entries are all unplaced (WGS family, ENV, TSA, GENE).
* ``empty_seq_id``, ``duplicate_seq_id``, ``seq_id_not_in_fasta``, ``no_entries``.
* ``unknown_category`` — not ``NONE``, ``GENE`` or a name in
  ``submission_categories.json``.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from ann_meta.model import Entry
from common.qualifiers import known_qualifier_names
from common.submission_category import get_category_rules

ROLES = ("chromosome", "plasmid", "organelle", "segment")
TYPES = ROLES + ("unplaced",)
STATUSES = ("complete", "partial")
TOPOLOGIES = ("circular", "linear", "")
# categories whose sequences may carry a role (chromosome / plasmid / ...);
# every other category takes unplaced entries only
ROLE_CATEGORIES = frozenset({"NONE", "GNM", "MAG", "TPA-GNM"})
RESERVED = frozenset({"ff_definition", "submitter_seqid", *ROLES})
_FORBIDDEN = re.compile(r"[\x00-\x1f\x7f]")


@dataclass
class Issue:
    field: str
    entry: str | None
    message: str


def _bad(value) -> bool:
    return isinstance(value, str) and bool(_FORBIDDEN.search(value))


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
    if category in ("NONE", "GENE"):
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
    for e in entries:
        allowed.update(q for q, _ in e.extra)

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
    return issues
