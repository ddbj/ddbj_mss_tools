"""Read the metadata of a DDBJ MSS .ann file into a :class:`MetaDoc`.

Rules (spec section 5):

* COMMON DBLINK / SUBMITTER / REFERENCE / COMMENT / DATE / ST_COMMENT go to
  ``common`` in the shape ``common.models.CommonModel`` accepts. Qualifiers
  without a field there go to ``common_extra[feature]``. Other COMMON
  features (DATATYPE, KEYWORD, DIVISION, unknown ones) are listed in
  ``untouched_common_features``.
* A COMMON ``source`` is expanded to every entry; an entry's own ``source``
  overrides it qualifier by qualifier.
* Per entry: ``chromosome`` / ``plasmid`` / ``organelle`` / ``segment`` give
  the role, ``submitter_seqid`` is dropped (it is regenerated), ``ff_definition``
  goes to ``Entry.ff_definition`` and a TOPOLOGY ``circular`` row to
  ``Entry.topology``. The remaining qualifiers go to ``source`` when every
  entry has the same value, otherwise to ``Entry.overrides`` (known names) or
  ``Entry.extra`` (unknown names).
* Anything that cannot be represented is reported as a :class:`Note`, never
  dropped silently. Empty values of known, non-flag qualifiers (DFAST's
  ``organism`` / ``strain`` placeholders) carry no information and are skipped.
"""
from __future__ import annotations

import re
from pathlib import Path

from ann_fix.annfile import AnnFeature, AnnFile
from ann_meta.category import EntrySketch, detect_category
from ann_meta.model import Entry, MetaDoc, Note
from common.models import ReferenceModel, SubmitterModel
from common.qualifiers import known_qualifier_names, qualifier_def

COMMON_TARGETS = ("DBLINK", "SUBMITTER", "REFERENCE", "COMMENT", "DATE", "ST_COMMENT")
ROLE_QUALIFIERS = ("chromosome", "plasmid", "organelle", "segment")
_DROPPED = ("submitter_seqid",)
_IDENTIFIER_EXCLUDED = {"organism", "entry", "submitter_seqid", *ROLE_QUALIFIERS}
_PLACEHOLDER = re.compile(r"@@\[([^\]]+)\]@@")
_COMPLETE = re.compile(r"\bcomplete\b", re.IGNORECASE)
_SIMPLE_LOCATION = re.compile(r"^1\.\.(\d+)$")

_DBLINK_NAMES = {
    "project": "project", "bioproject": "project",
    "biosample": "sample", "sample": "sample",
    "sequence read archive": "sequence_read_archive", "sequence_read_archive": "sequence_read_archive",
    "sra": "sequence_read_archive", "DRA": "sequence_read_archive",
}


# -- FASTA -------------------------------------------------------------------

def read_fasta_lengths(path: str | Path) -> dict[str, int]:
    """Return ``{id: length}``; the id is the first token after ``>``."""
    lengths: dict[str, int] = {}
    current: str | None = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line == "//":
                continue
            if line.startswith(">"):
                tokens = line[1:].split()
                current = tokens[0] if tokens else ""
                lengths[current] = 0
            elif current is not None:
                lengths[current] += len("".join(line.split()))
    return lengths


# -- COMMON ------------------------------------------------------------------

def _fold(quals: list[tuple[str, str]], allowed: set[str] | None, list_fields: set[str],
          rename: dict[str, str] | None = None) -> tuple[dict, list[tuple[str, str]]]:
    """Fold qualifier rows into a dict. Unknown names and repeats of scalar fields go to extra."""
    out: dict = {}
    extra: list[tuple[str, str]] = []
    for q, v in quals:
        key = rename.get(q, q) if rename else q
        if allowed is not None and key not in allowed:
            extra.append((q, v))
        elif key in list_fields:
            out.setdefault(key, []).append(v)
        elif key in out:
            extra.append((q, v))
        else:
            out[key] = v
    return out, extra


def _fold_dblink(quals: list[tuple[str, str]]) -> tuple[dict, list[tuple[str, str]]]:
    out, extra = _fold(quals, {"project", "sample", "sequence_read_archive"},
                       {"sample", "sequence_read_archive"}, _DBLINK_NAMES)
    if "sample" in out and len(out["sample"]) == 1:
        out["sample"] = out["sample"][0]
    return out, extra


def _read_common(feats: list[AnnFeature], rows, notes: list[Note]):
    common: dict = {}
    common_extra: dict[str, list[tuple[str, str]]] = {}
    untouched: list[str] = []
    raw: dict[str, list[tuple[str, str]]] = {}
    st_comments: list[dict] = []

    def add_extra(name: str, extra: list[tuple[str, str]]) -> None:
        if extra:
            common_extra.setdefault(name, []).extend(extra)

    for f in feats:
        quals = f.qualifiers(rows)
        if f.key == "source":
            continue
        raw.setdefault(f.key, []).extend(quals)
        if f.key not in COMMON_TARGETS:
            if f.key not in untouched:
                untouched.append(f.key)
            continue
        if f.key in ("DBLINK", "SUBMITTER", "DATE") and f.key in common:
            notes.append(Note("repeated_common_feature", None, f.key))
            add_extra(f.key, quals)
            continue
        if f.key == "DBLINK":
            value, extra = _fold_dblink(quals)
            common["DBLINK"] = value
        elif f.key == "SUBMITTER":
            value, extra = _fold(quals, set(SubmitterModel.model_fields), {"ab_name"})
            common["SUBMITTER"] = value
        elif f.key == "REFERENCE":
            value, extra = _fold(quals, set(ReferenceModel.model_fields), {"ab_name"})
            common.setdefault("REFERENCE", []).append(value)
        elif f.key == "COMMENT":
            lines = [v for q, v in quals if q == "line"]
            extra = [(q, v) for q, v in quals if q != "line"]
            value = {"line": lines[0] if len(lines) == 1 else lines} if lines else {}
            common.setdefault("COMMENT", []).append(value)
        elif f.key == "DATE":
            value, extra = _fold(quals, {"hold_date"}, set())
            common["DATE"] = value
        else:  # ST_COMMENT: any tag name is a field
            value, extra = _fold(quals, None, set())
            st_comments.append(value)
        add_extra(f.key, extra)
    if st_comments:
        common["ST_COMMENT"] = st_comments[0] if len(st_comments) == 1 else st_comments
    return common, common_extra, untouched, raw


# -- source ------------------------------------------------------------------

def _empty_placeholder(name: str, value: str) -> bool:
    d = qualifier_def(name)
    return value == "" and d is not None and not d.get("boolean")


def _qual_value(value: str) -> str | bool:
    return True if value == "" else value


def _note_unnamed(f: AnnFeature, rows, entry: str | None, notes: list[Note]) -> None:
    """Rows with a value but no qualifier name are invisible to AnnFeature.qualifiers."""
    for r in rows[f.first:f.last + 1]:
        if not r.cols[3] and r.cols[4]:
            notes.append(Note("unnamed_qualifier", entry, r.cols[4]))


def _topology(eid: str, feats: list[AnnFeature], rows, notes: list[Note]) -> str:
    """"circular" if a TOPOLOGY row says so; any other TOPOLOGY value gets a note."""
    topology = ""
    for f in feats:
        if f.key != "TOPOLOGY":
            continue
        quals = f.qualifiers(rows)
        if any("circular" in (q, v) for q, v in quals):
            topology = "circular"
        else:
            value = next((x for q, v in quals for x in (q, v) if x), "")
            notes.append(Note("topology_not_circular", eid, value))
    return topology


def extract(ann: AnnFile, fasta_lengths: dict[str, int] | None = None) -> MetaDoc:
    rows = ann.rows
    notes: list[Note] = []
    feats = ann.features()

    common_feats = [f for f in feats if f.entry == "COMMON"]
    common, common_extra, untouched, raw_common = _read_common(common_feats, rows, notes)
    for f in common_feats:
        if f.key in COMMON_TARGETS or f.key == "source":
            _note_unnamed(f, rows, None, notes)

    common_sources = [f for f in common_feats if f.key == "source"]
    common_src: list[tuple[str, str]] | None = None
    if common_sources:
        common_src = []
        for f in common_sources:
            common_src.extend(f.qualifiers(rows))
        if len(common_sources) > 1:
            notes.append(Note("repeated_source", None, ""))
        if common_sources[0].location != "1..E":
            notes.append(Note("source_location", None, common_sources[0].location))

    # entries in ann order
    entry_ids: list[str] = []
    by_entry: dict[str, list[AnnFeature]] = {}
    for f in feats:
        if f.entry and f.entry != "COMMON":
            if f.entry not in by_entry:
                entry_ids.append(f.entry)
                by_entry[f.entry] = []
            by_entry[f.entry].append(f)
            if f.key in ("source", "TOPOLOGY"):
                _note_unnamed(f, rows, f.entry, notes)
    entry_src: dict[str, list[tuple[str, str]]] = {}
    for eid in entry_ids:
        srcs = [f for f in by_entry[eid] if f.key == "source"]
        if not srcs:
            continue
        entry_src[eid] = [q for f in srcs for q in f.qualifiers(rows)]
        if len(srcs) > 1:
            notes.append(Note("repeated_source", eid, ""))
        length = (fasta_lengths or {}).get(eid)
        if length is not None and srcs[0].location != f"1..{length}":
            notes.append(Note("source_location", eid, srcs[0].location))

    has_common = common_src is not None
    layout = "common_source" if has_common and not entry_src else "per_entry"
    if has_common and entry_src:
        notes.append(Note("both_common_and_entry_source", None, ""))

    if fasta_lengths is not None:
        for eid in entry_ids:
            if eid not in fasta_lengths:
                notes.append(Note("entry_not_in_fasta", eid, ""))
        entry_ids += [i for i in fasta_lengths if i not in by_entry]

    if has_common and not entry_ids:
        notes.append(Note("common_source_without_entries", None, ""))

    if not has_common and not entry_src:
        notes.insert(0, Note("no_source", None, ""))
        category = detect_category(raw_common, [], None)
        return MetaDoc(category=category, common=common, common_extra=common_extra,
                       untouched_common_features=untouched, source={}, source_identifier=None,
                       entries=[], layout=layout, notes=notes)

    # merged qualifier rows per entry (COMMON source first, entry source overrides by name)
    merged: dict[str, list[tuple[str, str]]] = {}
    for eid in entry_ids:
        own = entry_src.get(eid)
        if own is None and not has_common:
            notes.append(Note("entry_without_source", eid, ""))
            continue
        if own is None:
            merged[eid] = list(common_src)
        elif has_common:
            names = {q for q, _ in own}
            merged[eid] = [(q, v) for q, v in common_src if q not in names] + own
        else:
            merged[eid] = list(own)

    known = known_qualifier_names()
    entries: list[Entry] = []
    sketches: list[EntrySketch] = []
    candidates: dict[str, dict[str, str | bool]] = {}   # eid -> single-valued qualifiers
    # eid -> ordered rows: (name, value, kind) kind in {"cand", "extra"}
    ordered: dict[str, list[tuple[str, str, str]]] = {}
    for eid in entry_ids:
        e = Entry(seq_id=eid, length=(fasta_lengths or {}).get(eid))
        entries.append(e)
        quals = merged.get(eid)
        has_topology = any(f.key == "TOPOLOGY" for f in by_entry.get(eid, []))
        e.topology = _topology(eid, by_entry.get(eid, []), rows, notes)
        if quals is None:
            sketches.append(EntrySketch(eid, has_topology, None, {}))
            continue
        counts: dict[str, int] = {}
        for q, _ in quals:
            counts[q] = counts.get(q, 0) + 1
        for q, c in counts.items():
            if c > 1 and q not in _DROPPED:
                notes.append(Note("repeated_qualifier", eid, q))
        cand: dict[str, str | bool] = {}
        seq: list[tuple[str, str, str]] = []
        role_taken = False
        for q, v in quals:
            if q in _DROPPED:
                continue
            if counts[q] > 1:
                seq.append((q, v, "extra"))
            elif q == "ff_definition":
                e.ff_definition = v
            elif q in ROLE_QUALIFIERS:
                if role_taken:
                    notes.append(Note("multiple_roles", eid, q))
                    seq.append((q, v, "extra"))
                else:
                    role_taken = True
                    e.type, e.seq_name = q, v
            elif _empty_placeholder(q, v):
                continue
            else:
                cand[q] = _qual_value(v)
                seq.append((q, v, "cand"))
        if e.ff_definition and _COMPLETE.search(e.ff_definition):
            e.status = "complete"
        candidates[eid] = cand
        ordered[eid] = seq
        sketch_quals: dict[str, str] = {}
        for q, v in quals:
            sketch_quals.setdefault(q, v)
        sketches.append(EntrySketch(eid, has_topology, e.ff_definition, sketch_quals))

    # qualifiers with the same value in every entry that has a source -> file-level source
    with_src = [eid for eid in entry_ids if eid in candidates]
    source: dict[str, str | bool] = {}
    for eid in with_src:
        for q, v in candidates[eid].items():
            if q in source:
                continue
            if all(q in candidates[o] and candidates[o][q] == v for o in with_src):
                source[q] = v
    for e in entries:
        if e.seq_id not in ordered:
            continue
        for q, v, kind in ordered[e.seq_id]:
            if kind == "extra":
                e.extra.append((q, _qual_value(v)))
            elif q in source:
                continue
            elif q in known:
                e.overrides[q] = _qual_value(v)
            else:
                e.extra.append((q, _qual_value(v)))

    source_identifier = None
    for e in entries:
        for x in _PLACEHOLDER.findall(e.ff_definition or ""):
            if x not in _IDENTIFIER_EXCLUDED:
                source_identifier = x
                break
        if source_identifier:
            break

    first_len = None
    if entries:
        first = entries[0]
        first_len = first.length
        if first_len is None:
            src = next((f for f in by_entry.get(first.seq_id, []) if f.key == "source"), None)
            m = _SIMPLE_LOCATION.match(src.location) if src else None
            first_len = int(m.group(1)) if m else None
    category = detect_category(raw_common, sketches, first_len)

    unknown = set(q for q in source if q not in known)
    for e in entries:
        unknown.update(q for q in e.overrides if q not in known)
        unknown.update(q for q, _ in e.extra if q not in known)

    return MetaDoc(category=category, common=common, common_extra=common_extra,
                   untouched_common_features=untouched, source=source,
                   source_identifier=source_identifier, entries=entries, layout=layout,
                   notes=notes, unknown_qualifiers=sorted(unknown))
