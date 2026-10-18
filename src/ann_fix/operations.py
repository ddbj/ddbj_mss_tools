"""Operations on a gff2mss .ann driven by ddbj-validator findings, plus the code -> default table.

Features are addressed by the identity of their first AnnRow object, never by line
number alone: rows are resolved from the untouched file before any edit, and a row
that an earlier operation removed or replaced simply no longer resolves (-> skip).
"""
from __future__ import annotations

import re

from .annfile import AnnFeature, AnnFile, AnnRow
from .report import Finding

ARTIFICIAL_LOCATION_DEFAULT = "low-quality sequence region"
ARTIFICIAL_LOCATION_VALUES = {"low-quality sequence region", "heterogeneous population sequenced"}

DEFAULTS: dict[str, tuple[str, dict]] = {
    "ANN5320": ("overlap-drop", {"keep": "cds"}),
    "ANN4300": ("misc", {}),
    "AXS6810": ("artificial-location", {"value": ARTIFICIAL_LOCATION_DEFAULT}),
    "ANN6840": ("artificial-location", {"value": ARTIFICIAL_LOCATION_DEFAULT}),   # intron < 10 bp
}
IGNORABLE: dict[str, str] = {
    "ANN2520": "同一遺伝子のアイソフォームが locus_tag を共有（想定どおり）",
    "ANN2542": "/gene 無しの locus_tag 共有（gff2mss の仕様。放置してよい）",
}
OPERATIONS = {"misc", "artificial-location", "drop", "overlap-drop", "overlap-misc"}
KEEP_VALUES = {"cds", "other"}

_DROP_ON_MISC = {"transl_table", "codon_start", "product", "gene", "inference", "transl_except",
                 "artificial_location"}
_PROTECTED_KEYS = {"source", "TOPOLOGY", "KEYWORD", "COMMON"}
_AL_KEYS = {"CDS", "mRNA"}
_TID = re.compile(r"submitter_transcript_id:\s*(.+?)\s*$")
_IN_CDS = re.compile(r"\bin CDS (\S+)")


# -- lookup helpers ----------------------------------------------------------

def _tids(ann: AnnFile, feat: AnnFeature) -> set[str]:
    for k, v in feat.qualifiers(ann.rows):
        if k == "note":
            m = _TID.search(v)
            if m:
                return {t.strip() for t in m.group(1).split(", ") if t.strip()}
    return set()


def _row(ann: AnnFile, feat: AnnFeature) -> AnnRow:
    return ann.rows[feat.first]


def _feat_of_row(ann: AnnFile, row: AnnRow) -> AnnFeature | None:
    for f in ann.features():
        if ann.rows[f.first] is row:
            return f
    return None


def _partners(ann: AnnFile, feat: AnnFeature) -> list[AnnFeature]:
    """mRNA(s) of a CDS / CDS of an mRNA, by submitter_transcript_id (fallback: locus_tag)."""
    if feat.key not in ("CDS", "mRNA"):
        return []
    want = "mRNA" if feat.key == "CDS" else "CDS"
    mine = _tids(ann, feat)
    cands = [f for f in ann.features() if f.entry == feat.entry and f.key == want]
    if mine:
        return [f for f in cands if _tids(ann, f) & mine]
    tag = feat.get(ann.rows, "locus_tag")
    return [f for f in cands if tag is not None and f.get(ann.rows, "locus_tag") == tag][:1]


def paired_mrna(ann: AnnFile, cds: AnnFeature) -> AnnFeature | None:
    p = _partners(ann, cds) if cds.key == "CDS" else []
    return p[0] if p else None


def paired_cds(ann: AnnFile, mrna: AnnFeature) -> AnnFeature | None:
    p = _partners(ann, mrna) if mrna.key == "mRNA" else []
    return p[0] if p else None


# -- operations --------------------------------------------------------------

def _delete_rows(ann: AnnFile, rows: list[AnnRow]) -> int:
    n = 0
    for r in rows:
        f = _feat_of_row(ann, r)
        if f is not None:
            ann.delete_feature(f)
            n += 1
    return n


def op_drop(ann: AnnFile, feat: AnnFeature, pair: bool = True, **kw) -> str:
    rows = [_row(ann, feat)] + ([_row(ann, p) for p in _partners(ann, feat)] if pair else [])
    n = _delete_rows(ann, rows)
    return f"dropped {feat.key} {feat.location}" + (f" and {n - 1} paired feature(s)" if n > 1 else "")


def op_misc(ann: AnnFile, feat: AnnFeature, reason: str = "", **kw) -> str:
    if feat.key == "mRNA":
        cds = paired_cds(ann, feat)
        if cds is not None:
            feat = cds
    target = _row(ann, feat)
    _delete_rows(ann, [_row(ann, p) for p in _partners(ann, feat)] if feat.key == "CDS" else [])
    feat = _feat_of_row(ann, target)
    quals = feat.qualifiers(ann.rows)
    tag = feat.get(ann.rows, "locus_tag")
    product = feat.get(ann.rows, "product")
    note = f"nonfunctional {feat.key} ({reason})" if reason else f"nonfunctional {feat.key}"
    if feat.key != "CDS" and product:
        note += f"; was {product}"
    new = ([("locus_tag", tag)] if tag is not None else []) + [("note", note)] + \
          [(k, v) for k, v in quals if k not in _DROP_ON_MISC and k != "locus_tag"]
    ann.replace_feature(feat, "misc_feature", feat.location, new)
    return f"{feat.key} {feat.location} -> misc_feature"


def op_artificial_location(ann: AnnFile, feat: AnnFeature, value: str = ARTIFICIAL_LOCATION_DEFAULT, **kw) -> str:
    rows = [_row(ann, feat)] + [_row(ann, p) for p in _partners(ann, feat)]
    added = 0
    for r in rows:
        t = _feat_of_row(ann, r)
        if t is not None and t.key in _AL_KEYS and ("artificial_location", value) not in t.qualifiers(ann.rows):
            ann.add_qualifier(t, "artificial_location", value)
            added += 1
    return f"artificial_location on {feat.key} {feat.location} ({added} added)"


def _other_side(ann: AnnFile, feat: AnnFeature, finding: Finding | None) -> AnnFeature | None:
    """For an overlap finding on *feat*, the feature named in the message ("in CDS <loc>")."""
    if finding is None:
        return None
    m = _IN_CDS.search(finding.message)
    if not m:
        return None
    loc = m.group(1)
    while loc.endswith(")") and loc.count(")") > loc.count("("):
        loc = loc[:-1]
    return next((f for f in ann.features() if f.entry == feat.entry and f.key == "CDS" and f.location == loc), None)


def _overlap_loser(ann, feat, finding, keep):
    if keep == "cds":
        return feat
    other = _other_side(ann, feat, finding)
    if other is None:
        raise ValueError(f"cannot locate the CDS named in the finding at line {finding.line if finding else '?'}")
    return other


def op_overlap_drop(ann: AnnFile, feat: AnnFeature, finding: Finding | None = None, keep: str = "cds", **kw) -> str:
    return op_drop(ann, _overlap_loser(ann, feat, finding, keep))


def op_overlap_misc(ann: AnnFile, feat: AnnFeature, finding: Finding | None = None, keep: str = "cds", **kw) -> str:
    return op_misc(ann, _overlap_loser(ann, feat, finding, keep), reason=finding.code if finding else "overlap")


_OPS = {"misc": op_misc, "artificial-location": op_artificial_location, "drop": op_drop,
        "overlap-drop": op_overlap_drop, "overlap-misc": op_overlap_misc}


def _resolve(ann: AnnFile, f: Finding) -> AnnRow | None:
    """The first row of the feature a finding points at, or None when the report does not fit this file."""
    if f.line is None or not 1 <= f.line <= len(ann.rows):
        return None
    row = ann.rows[f.line - 1]
    if row.cols[1] != f.feature:
        return None
    if "Found:" in f.message and row.cols[2] and row.cols[2] not in f.message:
        return None
    return row


def apply(ann: AnnFile, findings: list[Finding], choice: dict[str, tuple[str, dict]]) -> list[str]:
    log: list[str] = []
    present = {f.code for f in findings}
    for code in choice:
        if code not in present:
            log.append(f"no findings for {code}")
    todo = []
    unresolved = 0
    for f in findings:
        if f.code not in choice:
            continue
        if f.line is None:
            continue                      # file-level findings carry no feature
        row = _resolve(ann, f)
        if row is None:
            unresolved += 1
            continue
        todo.append((f, row))
    skipped_file_level = {f.code for f in findings if f.code in choice and f.line is None}
    for code in sorted(skipped_file_level):
        log.append(f"skip: {code} is file-level (no feature)")
    if unresolved:
        raise ValueError(f"report does not match this .ann ({unresolved} finding(s) do not point at the "
                         f"reported feature); re-run ddbj-validator on this file first")
    for f, row in todo:
        feat = _feat_of_row(ann, row)
        if feat is None:
            log.append(f"skip: line {f.line} {f.feature} already changed by an earlier operation ({f.code})")
            continue
        op, kwargs = choice[f.code]
        if op == "artificial-location" and feat.key not in _AL_KEYS:
            log.append(f"skip: line {f.line} {feat.key}: artificial_location is only for CDS/mRNA ({f.code})")
            continue
        if op != "artificial-location" and feat.key in _PROTECTED_KEYS:
            log.append(f"skip: line {f.line} {feat.key}: {op} is not allowed on this feature ({f.code})")
            continue
        reason = f"{f.code}: {f.message.split('(')[0].strip()}"
        try:
            msg = _OPS[op](ann, feat, finding=f, reason=reason, **kwargs)
        except ValueError as exc:
            log.append(f"skip: line {f.line} {feat.key}: {exc} ({f.code})")
            continue
        log.append(f"{f.code} line {f.line}: {msg}")
    return log
