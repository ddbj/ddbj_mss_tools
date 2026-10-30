"""Write edited metadata back into a DDBJ MSS .ann file (`ann_meta apply`).

Only these rows are rewritten; every other row keeps its bytes (columns and
line ending), so annotation features survive unchanged:

1. COMMON DBLINK / SUBMITTER / REFERENCE / COMMENT / DATE / ST_COMMENT. The
   existing ones are removed and the rows ``create_common`` builds from the
   common JSON are inserted where the first of them was (or at the end of
   COMMON when there was none), in that canonical order. Qualifiers under
   ``common["_extra"][feature]`` follow the generated rows of that feature.
   DATATYPE, DIVISION, KEYWORD and unknown COMMON features are not touched,
   even when the common JSON has keys of those names.
2. ``source`` and ``TOPOLOGY``. Layout ``common_source`` (one COMMON source at
   ``1..E``) is used only when the category is WGS-like and every entry is
   unplaced without per-entry qualifiers; otherwise every entry gets its own
   source (layout ``per_entry``). In ``common_source`` a contig whose only rows
   were its source has no rows afterwards (as in mss_builder's WGS output: the
   COMMON source covers every FASTA entry). A per-entry source is ``SOURCE`` -> overrides -> role
   qualifier -> ``ff_definition`` -> ``extra``; a ``TOPOLOGY circular`` row
   goes directly before it when the entry is circular. Existing TOPOLOGY rows
   and extra source features of the entry are removed.

The entry table must cover every entry of the .ann: an entry missing from it
raises ``ValueError("entry_missing_from_table: <id>")`` before anything is
changed. ``apply`` does not otherwise validate; run
:func:`ann_meta.validate.validate` first.
"""
from __future__ import annotations

import dataclasses

from ann_fix.annfile import AnnFeature, AnnFile
from ann_meta.extract import COMMON_TARGETS
from ann_meta.model import Entry
from common.common_builder import create_common
from common.qualifiers import qualifier_def
from common.source_builder import (
    _FLAG_QUALIFIERS,
    SequenceRoleEntry,
    _flag_is_set,
    create_source_feature,
    ff_definition,
    source_qualifier,
)
from common.submission_category import SubmissionCategoryRules, get_category_rules

# Categories whose all-unplaced entries share one COMMON source (WGS, MAG-WGS = required_datatype WGS).
# ENV has no required_datatype in submission_categories.json; it is listed by the spec (contig sets like WGS).
COMMON_SOURCE_CATEGORIES = frozenset({"WGS", "MAG-WGS", "ENV"})
# Categories whose ff_definition is regenerated from the role; for the others
# (NONE, GENE) the entry's ff_definition is written as given.
REGENERATE_FF_CATEGORIES = frozenset({"GNM", "WGS", "MAG", "MAG-WGS", "ENV", "TSA", "TSA-ENV"})
_NO_RULES = ("NONE", "GENE")

Quals = list[tuple[str, str]]


def _rules(category: str) -> SubmissionCategoryRules:
    return SubmissionCategoryRules() if category in _NO_RULES else get_category_rules(category)


def _is_flag(name: str) -> bool:
    d = qualifier_def(name)
    return name in _FLAG_QUALIFIERS or bool(d and d.get("boolean"))


def _qual_rows(name: str, value) -> Quals:
    """(name, value) -> qualifier rows. Flags: on -> valueless row, off -> nothing."""
    if value is None or value is False:
        return []
    if _is_flag(name):
        return [(name, "")] if _flag_is_set(value) else []
    if value is True:
        return [(name, "")]
    if isinstance(value, list):
        return [(name, str(v)) for v in value]
    return [(name, str(value))]


def _blocks(rows: list[list[str]]) -> list[tuple[str, str, Quals]]:
    """Split builder rows (5-element lists) into (key, location, qualifiers) per feature."""
    out: list[tuple[str, str, Quals]] = []
    for r in rows:
        if r[1]:
            out.append((r[1], r[2], []))
        if r[3] or r[4]:
            out[-1][2].append((r[3], r[4]))
    return out


def _common_end(ann: AnnFile) -> int:
    common = [f for f in ann.features() if f.entry == "COMMON"]
    return common[-1].last + 1 if common else 0


def _common_marker(ann: AnnFile, at: int) -> str:
    """'COMMON' when the block inserted at *at* must carry the COMMON marker."""
    return "" if any(r.cols[0] == "COMMON" for r in ann.rows[:at]) else "COMMON"


def _insert_blocks(ann: AnnFile, at: int, blocks: list[tuple[str, str, Quals]]) -> None:
    for key, location, quals in blocks:
        f = ann.insert_feature(at, _common_marker(ann, at), key, location, quals)
        at = f.last + 1


# -- COMMON --------------------------------------------------------------------

def _apply_common_targets(ann: AnnFile, common: dict) -> None:
    targets = [f for f in ann.features() if f.entry == "COMMON" and f.key in COMMON_TARGETS]
    at = targets[0].first if targets else None
    for f in reversed(targets):
        ann.delete_feature(f)
    if at is None:
        at = _common_end(ann)
    extra = common.get("_extra") or {}
    blocks: list[tuple[str, str, Quals]] = []
    for name in COMMON_TARGETS:
        own = _blocks(create_common({name: common[name]})) if common.get(name) else []
        pairs = [(str(k), "" if v is None or v is True else str(v)) for k, v in extra.get(name, [])]
        if pairs:
            if own:
                own[0][2].extend(pairs)
            else:
                own = [(name, "", pairs)]
        blocks.extend(own)
    _insert_blocks(ann, at, blocks)


def _common_source_blocks(common: dict, category: str, source_identifier: str) -> list[tuple[str, str, Quals]]:
    source = {}
    for k, v in (common.get("SOURCE") or {}).items():
        if _is_flag(k):
            source[k] = v
        elif v is True:
            source[k] = ""
        elif v is not False and v is not None:
            source[k] = v
    rows = create_source_feature(category, None, None, None, source,
                                 source_modifier_key=source_identifier, use_meta_expression=True)
    return _blocks(rows)


# -- entries -------------------------------------------------------------------

def _decide_layout(category: str, entries: list[Entry]) -> str:
    if category not in COMMON_SOURCE_CATEGORIES:
        return "per_entry"
    if any(e.type != "unplaced" or e.overrides or e.extra for e in entries):
        return "per_entry"
    return "common_source"


class _SourceBuilder:
    def __init__(self, common: dict, entries: list[Entry], category: str, source_identifier: str):
        self.base = dict(common.get("SOURCE") or {})
        self.rules = _rules(category)
        self.regenerate = category in REGENERATE_FF_CATEGORIES
        self.is_wgs = category in COMMON_SOURCE_CATEGORIES and all(e.type == "unplaced" for e in entries)
        self.source_identifier = source_identifier or None
        self.chromosomes = sum(1 for e in entries if e.type == "chromosome")
        self.segments = sum(1 for e in entries if e.type == "segment")

    def qualifiers(self, e: Entry) -> Quals:
        merged = dict(self.base)
        for k in self.rules.auto_source_qualifiers:
            merged.setdefault(k, True)
        merged.update(e.overrides)
        rows = [r for k, v in merged.items() for r in _qual_rows(k, v)]
        role = SequenceRoleEntry(e.seq_id, e.type, e.seq_name, e.status, e.topology == "circular")
        for k, v in source_qualifier(role, e.seq_id, self.is_wgs, self.segments).items():
            rows.extend(_qual_rows(k, v))
        if self.regenerate:
            mol_type = merged.get("mol_type")
            ff = ff_definition(role, self.source_identifier,
                               mol_type if isinstance(mol_type, str) else "genomic DNA",
                               self.is_wgs, self.chromosomes, self.segments)
        else:
            ff = e.ff_definition or ""
        if ff:
            rows.append(("ff_definition", ff))
        for k, v in e.extra:
            rows.extend(_qual_rows(k, v))
        return rows


def _rewrite_entry(ann: AnnFile, e: Entry, feats: list[AnnFeature], layout: str, build: _SourceBuilder) -> None:
    srcs = [f for f in feats if f.key == "source"]
    keep = srcs[0] if layout == "per_entry" and srcs else None
    doomed = [f for f in feats if f.key == "TOPOLOGY" or (f.key == "source" and f is not keep)]
    entry_first = feats[0].first
    for f in sorted(doomed, key=lambda f: f.first, reverse=True):
        ann.delete_feature(f)
        if keep is not None and f.first < keep.first:
            n = f.last - f.first + 1
            keep = dataclasses.replace(keep, first=keep.first - n, last=keep.last - n)
    emptied = len(doomed) == len(feats)

    def marker(at: int) -> str:
        return e.seq_id if emptied or (at < len(ann.rows) and ann.rows[at].cols[0] == e.seq_id) else ""

    if layout == "per_entry":
        quals = build.qualifiers(e)
        if keep is not None:
            ann.replace_feature(keep, "source", keep.location, quals)
            at = keep.first
        else:
            at = entry_first
            ann.insert_feature(at, marker(at), "source", _location(e), quals)
            emptied = False
    else:
        at = entry_first
    if e.topology == "circular":
        ann.insert_feature(at, marker(at), "TOPOLOGY", "", [("circular", "")])


def _location(e: Entry) -> str:
    if e.length is None:
        raise ValueError(f"cannot write a source for {e.seq_id!r}: sequence length unknown")
    return f"1..{e.length}"


def _append_entry(ann: AnnFile, e: Entry, build: _SourceBuilder) -> None:
    if e.topology == "circular":
        ann.insert_feature(len(ann.rows), e.seq_id, "TOPOLOGY", "", [("circular", "")])
        ann.insert_feature(len(ann.rows), "", "source", _location(e), build.qualifiers(e))
    else:
        ann.insert_feature(len(ann.rows), e.seq_id, "source", _location(e), build.qualifiers(e))


# -- entry point -----------------------------------------------------------------

def apply(ann: AnnFile, common: dict, entries: list[Entry], *, category: str) -> AnnFile:
    """Rewrite *ann* in place from *common* (common JSON) and *entries*; return it."""
    rules = _rules(category)
    source_identifier = (common.get("SOURCE_IDENTIFIER") or common.get("INFRASPECIFIC_NAME_MODIFIER")
                         or rules.source_identifier or "")

    listed = {e.seq_id for e in entries}
    for f in ann.features():
        if f.entry and f.entry != "COMMON" and f.entry not in listed:
            raise ValueError(f"entry_missing_from_table: {f.entry}")

    _apply_common_targets(ann, common)

    layout = _decide_layout(category, entries)

    for f in reversed([f for f in ann.features() if f.entry == "COMMON" and f.key == "source"]):
        ann.delete_feature(f)
    if layout == "common_source":
        _insert_blocks(ann, _common_end(ann), _common_source_blocks(common, category, source_identifier))

    by_entry: dict[str, list[AnnFeature]] = {}
    for f in ann.features():
        if f.entry and f.entry != "COMMON":
            by_entry.setdefault(f.entry, []).append(f)
    build = _SourceBuilder(common, entries, category, source_identifier)
    present = [e for e in entries if e.seq_id in by_entry]
    for e in sorted(present, key=lambda e: by_entry[e.seq_id][0].first, reverse=True):
        _rewrite_entry(ann, e, by_entry[e.seq_id], layout, build)
    if layout == "per_entry":           # in common_source the COMMON source covers FASTA-only entries
        for e in entries:
            if e.seq_id not in by_entry:
                _append_entry(ann, e, build)
    return ann
