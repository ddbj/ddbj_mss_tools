"""Parse DDBJ MSS annotation files (.ann or .annt.tsv) into structured data."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class DbLink:
    project: str = ""
    biosample: str = ""
    sra: list[str] = field(default_factory=list)   # Sequence Read Archive accessions
    extra: list[tuple[str, str]] = field(default_factory=list)  # other DBLINK entries


@dataclass
class Submitter:
    ab_names: list[str] = field(default_factory=list)
    consrtm: str = ""
    contact: str = ""
    email: str = ""
    url: str = ""
    institute: str = ""
    department: str = ""
    country: str = ""
    state: str = ""
    city: str = ""
    street: str = ""
    zip: str = ""


@dataclass
class Reference:
    title: str = ""
    ab_names: list[str] = field(default_factory=list)
    consrtm: str = ""
    status: str = "Unpublished"
    year: str = ""
    journal: str = ""
    volume: str = ""
    from_page: str = ""
    to_page: str = ""
    doi: str = ""
    pubmed: str = ""
    remark: str = ""


@dataclass
class StComment:
    tagset_id: str = ""
    fields: list[tuple[str, str]] = field(default_factory=list)


@dataclass
class CommonBlock:
    dblink: DbLink = field(default_factory=DbLink)
    submitter: Submitter = field(default_factory=Submitter)
    references: list[Reference] = field(default_factory=list)
    comment_blocks: list[list[str]] = field(default_factory=list)  # each inner list = one COMMENT block
    st_comments: list[StComment] = field(default_factory=list)
    hold_date: str = ""
    keywords: list[str] = field(default_factory=list)
    source_feature: Optional[Feature] = None  # source feature defined in COMMON block


@dataclass
class Feature:
    name: str
    location: str
    qualifiers: list[tuple[str, str]] = field(default_factory=list)

    def get_qualifier(self, key: str) -> Optional[str]:
        for k, v in self.qualifiers:
            if k == key:
                return v
        return None

    def get_qualifiers(self, key: str) -> list[str]:
        return [v for k, v in self.qualifiers if k == key]

    def has_qualifier(self, key: str) -> bool:
        return any(k == key for k, _ in self.qualifiers)


@dataclass
class EntryBlock:
    entry_id: str
    topology: str = "linear"
    features: list[Feature] = field(default_factory=list)
    keywords: list[str] = field(default_factory=list)

    @property
    def source_feature(self) -> Optional[Feature]:
        return next((f for f in self.features if f.name == "source"), None)


def _pad(row: list[str], n: int = 5) -> list[str]:
    return (row + [""] * n)[:n]


def parse_ann(filepath: str | Path) -> tuple[CommonBlock, list[EntryBlock]]:
    """Parse a DDBJ MSS annotation file. Returns (common, entries)."""
    with open(filepath, newline="", encoding="utf-8") as fh:
        rows = [_pad(r) for r in csv.reader(fh, delimiter="\t")]

    common = CommonBlock()
    entries: list[EntryBlock] = []

    i = 0
    n = len(rows)

    while i < n:
        entry_col = rows[i][0]
        if entry_col == "COMMON":
            i = _parse_common(rows, i, common)
        elif entry_col:
            entry = EntryBlock(entry_id=entry_col)
            i = _parse_entry(rows, i, entry)
            entries.append(entry)
        else:
            i += 1

    return common, entries


# ── COMMON block ──────────────────────────────────────────────────────────────

def _parse_common(rows: list[list[str]], i: int, common: CommonBlock) -> int:
    """Parse COMMON block starting at row i. Returns next row index."""
    n = len(rows)
    # row i has entry_col=="COMMON"; col[1] is the first section name
    current_section = rows[i][1]
    current_ref: Optional[Reference] = None
    current_st: Optional[StComment] = None
    current_comment: Optional[list[str]] = None

    def flush_section():
        nonlocal current_ref, current_st, current_comment
        if current_ref is not None:
            common.references.append(current_ref)
            current_ref = None
        if current_st is not None:
            common.st_comments.append(current_st)
            current_st = None
        if current_comment is not None:
            common.comment_blocks.append(current_comment)
            current_comment = None

    def start_section(name: str, loc: str = ""):
        nonlocal current_section, current_ref, current_st, current_comment
        flush_section()
        current_section = name
        if name == "REFERENCE":
            current_ref = Reference()
        elif name == "ST_COMMENT":
            current_st = StComment()
        elif name == "COMMENT":
            current_comment = []
        elif name == "source":
            common.source_feature = Feature(name="source", location=loc)

    def process_row(section: str, key: str, val: str):
        nonlocal current_ref, current_st, current_comment
        if section == "DBLINK":
            k = key.lower()
            if k == "project":
                common.dblink.project = val
            elif k == "biosample":
                common.dblink.biosample = val
            elif k in ("sequence read archive", "dra", "sra"):
                common.dblink.sra.append(val)
            else:
                common.dblink.extra.append((key, val))
        elif section == "SUBMITTER":
            if key == "ab_name":
                common.submitter.ab_names.append(val)
            elif key == "consrtm":
                common.submitter.consrtm = val
            else:
                setattr(common.submitter, key.replace(" ", "_"), val)
        elif section == "REFERENCE" and current_ref is not None:
            if key == "ab_name":
                current_ref.ab_names.append(val)
            elif key == "consrtm":
                current_ref.consrtm = val
            elif key == "title":
                current_ref.title = val
            elif key == "status":
                current_ref.status = val
            elif key == "year":
                current_ref.year = val
            elif key == "journal":
                current_ref.journal = val
            elif key == "volume":
                current_ref.volume = val
            elif key in ("from_page", "start_page"):
                current_ref.from_page = val
            elif key in ("to_page", "end_page"):
                current_ref.to_page = val
            elif key == "doi":
                current_ref.doi = val
            elif key == "pubmed":
                current_ref.pubmed = val
        elif section == "COMMENT" and current_comment is not None:
            if key == "line":
                current_comment.append(val)
        elif section == "ST_COMMENT" and current_st is not None:
            if key == "tagset_id":
                current_st.tagset_id = val
            else:
                current_st.fields.append((key, val))
        elif section == "source" and common.source_feature is not None and key:
            common.source_feature.qualifiers.append((key, val))
        elif section == "DATE":
            if key == "hold_date":
                common.hold_date = val
        elif section == "KEYWORD":
            common.keywords.append(val)

    # First row
    process_row(current_section, rows[i][3], rows[i][4])
    i += 1

    while i < n:
        entry_col, feat_col, loc_col, key_col, val_col = rows[i]
        if entry_col:  # new top-level block
            break
        if feat_col:
            start_section(feat_col, loc_col)
        process_row(current_section, key_col, val_col)
        i += 1

    flush_section()
    return i


# ── Entry block ───────────────────────────────────────────────────────────────

def _parse_entry(rows: list[list[str]], i: int, entry: EntryBlock) -> int:
    """Parse one entry block starting at row i. Returns next row index."""
    n = len(rows)
    current_feat: Optional[Feature] = None

    def start_feature(feat_name: str, location: str, key: str, val: str):
        nonlocal current_feat
        if feat_name == "TOPOLOGY":
            entry.topology = "circular" if key == "circular" else "linear"
            current_feat = None
        elif feat_name == "KEYWORD":
            entry.keywords.append(key or val)
            current_feat = None
        else:
            current_feat = Feature(name=feat_name, location=location)
            if key:
                current_feat.qualifiers.append((key, val))
            entry.features.append(current_feat)

    # First row
    start_feature(rows[i][1], rows[i][2], rows[i][3], rows[i][4])
    i += 1

    while i < n:
        entry_col, feat_col, loc_col, key_col, val_col = rows[i]
        if entry_col:  # new entry
            break
        if feat_col:
            start_feature(feat_col, loc_col, key_col, val_col)
        elif current_feat is not None and key_col:
            current_feat.qualifiers.append((key_col, val_col))
        i += 1

    return i
