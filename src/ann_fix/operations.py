"""Operations on a gff2mss .ann driven by ddbj-validator findings, plus the code -> default table."""
from __future__ import annotations

import re

from .annfile import AnnFeature, AnnFile
from .report import Finding

ARTIFICIAL_LOCATION_DEFAULT = "low-quality sequence region"

DEFAULTS: dict[str, tuple[str, dict]] = {
    "ANN5320": ("overlap-drop", {"keep": "cds"}),
    "ANN4300": ("misc", {}),
    "AXS6810": ("artificial-location", {"value": ARTIFICIAL_LOCATION_DEFAULT}),
}
IGNORABLE: dict[str, str] = {
    "ANN2520": "同一遺伝子のアイソフォームが locus_tag を共有（想定どおり）",
    "ANN2542": "/gene 無しの locus_tag 共有（gff2mss の仕様。放置してよい）",
}
OPERATIONS = {"misc", "artificial-location", "drop", "overlap-drop", "overlap-misc"}

_DROP_ON_MISC = {"transl_table", "codon_start", "product", "gene", "inference", "transl_except",
                 "artificial_location"}
_TID = re.compile(r"submitter_transcript_id:\s*(.+?)\s*$")
_IN_CDS = re.compile(r"\bin CDS (\S+)")


# -- lookup helpers ----------------------------------------------------------

def _tid(ann: AnnFile, feat: AnnFeature) -> str | None:
    for k, v in feat.qualifiers(ann.rows):
        if k == "note":
            m = _TID.search(v)
            if m:
                return m.group(1).rstrip(",")
    return None


def _ident(ann: AnnFile, feat: AnnFeature) -> tuple:
    return (feat.entry, feat.key, feat.location, feat.get(ann.rows, "locus_tag"))


def _find(ann: AnnFile, ident: tuple) -> AnnFeature | None:
    """Locate a feature by identity captured before rows moved."""
    entry, key, loc, tag = ident
    for f in ann.features():
        if (f.entry, f.key, f.location) == (entry, key, loc) and f.get(ann.rows, "locus_tag") == tag:
            return f
    return None


def _partner(ann: AnnFile, feat: AnnFeature, key: str) -> AnnFeature | None:
    tid = _tid(ann, feat)
    tag = feat.get(ann.rows, "locus_tag")
    cands = [f for f in ann.features() if f.entry == feat.entry and f.key == key]
    if tid is not None:
        for f in cands:
            if _tid(ann, f) == tid:
                return f
        return None
    for f in cands:
        if tag is not None and f.get(ann.rows, "locus_tag") == tag:
            return f
    return None


def paired_mrna(ann: AnnFile, cds: AnnFeature) -> AnnFeature | None:
    return _partner(ann, cds, "mRNA") if cds.key == "CDS" else None


def paired_cds(ann: AnnFile, mrna: AnnFeature) -> AnnFeature | None:
    return _partner(ann, mrna, "CDS") if mrna.key == "mRNA" else None


def _pair(ann: AnnFile, feat: AnnFeature) -> AnnFeature | None:
    return paired_mrna(ann, feat) if feat.key == "CDS" else paired_cds(ann, feat)


# -- operations --------------------------------------------------------------

def op_drop(ann: AnnFile, feat: AnnFeature, pair: bool = True, **kw) -> str:
    partner = _pair(ann, feat) if pair else None
    pid = _ident(ann, partner) if partner is not None else None
    ann.delete_feature(feat)
    if pid is not None:
        partner = _find(ann, pid)
        if partner is not None:
            ann.delete_feature(partner)
    return f"dropped {feat.key} {feat.location}" + (" and its pair" if partner else "")


def op_misc(ann: AnnFile, feat: AnnFeature, reason: str = "", **kw) -> str:
    if feat.key == "mRNA":
        cds = paired_cds(ann, feat)
        if cds is not None:
            feat = cds
    partner = _pair(ann, feat) if feat.key == "CDS" else None
    if partner is not None:
        fid = _ident(ann, feat)
        ann.delete_feature(partner)
        feat = _find(ann, fid)
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
    targets = [_ident(ann, feat)]
    partner = _pair(ann, feat)
    if partner is not None:
        targets.append(_ident(ann, partner))
    added = 0
    for ident in targets:
        t = _find(ann, ident)
        if t is not None and ("artificial_location", value) not in t.qualifiers(ann.rows):
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


def apply(ann: AnnFile, findings: list[Finding], choice: dict[str, tuple[str, dict]]) -> list[str]:
    log: list[str] = []
    present = {f.code for f in findings}
    for code in choice:
        if code not in present:
            log.append(f"no findings for {code}")
    todo = [f for f in findings if f.code in choice]
    # ponytail: descending line order keeps earlier indices valid; a pair deleted above a
    # not-yet-processed feature shifts it and that finding is then skipped (logged), not misapplied.
    todo.sort(key=lambda f: -(f.line or 0))
    for f in todo:
        if f.line is None:
            log.append(f"skip: {f.code} is file-level (no line)")
            continue
        feat = ann.feature_at_line(f.line)
        if feat is None or feat.key != f.feature:
            log.append(f"skip: line {f.line} no longer matches {f.feature} ({f.code})")
            continue
        op, kwargs = choice[f.code]
        reason = f"{f.code}: {f.message.split('(')[0].strip()}"
        msg = _OPS[op](ann, feat, finding=f, reason=reason, **kwargs)
        log.append(f"{f.code} line {f.line}: {msg}")
    return log
