"""ann_meta: read (extract) and rewrite (apply) the metadata of a DDBJ MSS .ann.

``ann_meta extract --ann X [--fasta Y] [--json] [--entries-out E.tsv]``
    Print ``MetaDoc.to_json()`` on stdout (always JSON; ``--json`` is accepted
    for symmetry) and optionally write the entry table TSV. Exit 0, or 1 when
    the .ann / FASTA cannot be read.

``ann_meta apply --ann X --common C.json --entries E.tsv --out Z --category CAT
[--fasta Y] [--allow-unknown q1,q2]``
    Check the inputs; on any problem print ``{"issues": [{"field", "entry",
    "message"}]}`` on stdout and exit 2 without writing. Otherwise rewrite the
    COMMON / source rows and write Z through a temporary file in Z's directory
    (so ``--out`` may equal ``--ann``). Exit 0. Unreadable inputs exit 1.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tempfile
from pathlib import Path

from pydantic import ValidationError

from ann_fix.annfile import AnnFile
from ann_meta.apply import apply
from ann_meta.entries_tsv import read_entries, write_entries
from ann_meta.extract import COMMON_TARGETS, extract, read_fasta_lengths
from ann_meta.validate import Issue, validate
from common.models import CommonModel

_RANGE = re.compile(r"^1\.\.(\d+)$")


class _UserError(Exception):
    """Unreadable input: exit 1 with a one-line message."""


def _err(msg: str) -> None:
    sys.stderr.write(f"error: {msg}\n")


def _read_ann(path: str) -> AnnFile:
    try:
        return AnnFile.read(path)
    except (OSError, UnicodeDecodeError) as exc:
        raise _UserError(f"cannot read .ann {path}: {exc}") from None


def _read_fasta(path: str) -> dict[str, int]:
    try:
        return read_fasta_lengths(path)
    except OSError as exc:
        raise _UserError(f"cannot read FASTA {path}: {exc}") from None


def _issues_out(issues: list[Issue]) -> int:
    print(json.dumps({"issues": [{"field": i.field, "entry": i.entry, "message": i.message} for i in issues]},
                     ensure_ascii=False, indent=2))
    return 2


# -- extract ---------------------------------------------------------------------

def _cmd_extract(args) -> int:
    ann = _read_ann(args.ann)
    lengths = _read_fasta(args.fasta) if args.fasta else None
    doc = extract(ann, lengths)
    if args.entries_out:
        try:
            write_entries(doc.entries, args.entries_out)
        except OSError as exc:
            raise _UserError(f"cannot write {args.entries_out}: {exc}") from None
    print(json.dumps(doc.to_json(), ensure_ascii=False, indent=2))
    return 0


# -- apply -----------------------------------------------------------------------

def _load_common(path: str) -> dict:
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        raise _UserError(f"cannot read common JSON {path}: {exc}") from None
    if not isinstance(data, dict):
        raise _UserError(f"common JSON {path}: the top level must be an object")
    return data


def _shape_issues(common: dict) -> list[Issue]:
    issues: list[Issue] = []
    try:
        CommonModel.model_validate(common)
    except ValidationError as exc:
        for e in exc.errors():
            field = ".".join(str(p) for p in e["loc"]) or "common"
            issues.append(Issue(field, None, f"bad_shape: {e['msg']}"))
    extra = common.get("_extra")
    if extra is None:
        return issues
    if not isinstance(extra, dict):
        return issues + [Issue("_extra", None, "bad_shape: must be an object of feature -> [[name, value], ...]")]
    for feature, pairs in extra.items():
        if feature not in COMMON_TARGETS:   # apply would drop it silently
            issues.append(Issue(f"_extra.{feature}", None, "bad_shape: not a COMMON target"))
            continue
        if not isinstance(pairs, list):
            issues.append(Issue(f"_extra.{feature}", None, "bad_shape: must be a list of [name, value] pairs"))
            continue
        for i, pair in enumerate(pairs):
            ok = (isinstance(pair, list) and len(pair) == 2 and isinstance(pair[0], str)
                  and isinstance(pair[1], (str, bool)))
            if not ok:
                issues.append(Issue(f"_extra.{feature}.{i}", None,
                                    "bad_shape: must be [name, value] with a string name and a string/boolean value"))
    return issues


def _source_lengths(ann: AnnFile) -> dict[str, int]:
    """Entry -> N for entries whose source is at ``1..N``."""
    out: dict[str, int] = {}
    for f in ann.features():
        if f.key == "source" and f.entry and f.entry != "COMMON":
            m = _RANGE.match(f.location.strip())
            if m:
                out.setdefault(f.entry, int(m.group(1)))
    return out


def _write_atomic(ann: AnnFile, out: str) -> None:
    target = Path(out)
    fd, tmp = tempfile.mkstemp(dir=target.parent,
                               prefix=f".{target.name}.", suffix=".tmp")
    os.close(fd)
    try:
        ann.write(tmp)
        os.replace(tmp, target)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _cmd_apply(args) -> int:
    common = _load_common(args.common)
    shape = _shape_issues(common)
    if shape:
        return _issues_out(shape)
    try:
        entries = read_entries(args.entries)
    except (OSError, UnicodeDecodeError) as exc:
        raise _UserError(f"cannot read entry table {args.entries}: {exc}") from None
    except ValueError as exc:
        return _issues_out([Issue("entries", None, f"bad_entries_tsv: {exc}")])
    ann = _read_ann(args.ann)

    if args.fasta:
        lengths = _read_fasta(args.fasta)
        fasta_ids: list[str] | None = list(lengths)
    else:
        lengths = _source_lengths(ann)
        fasta_ids = None
    for e in entries:
        e.length = lengths.get(e.seq_id)

    if args.allow_unknown is None:
        allowed = set(extract(ann).unknown_qualifiers)
    else:
        allowed = {q.strip() for q in args.allow_unknown.split(",") if q.strip()}

    issues = validate(common, entries, category=args.category, fasta_ids=fasta_ids, allowed_unknown=allowed)
    if issues:
        return _issues_out(issues)
    try:
        apply(ann, common, entries, category=args.category)
    except ValueError as exc:
        msg = str(exc)
        if msg.startswith("entry_missing_from_table:"):
            return _issues_out([Issue("entries", msg.split(":", 1)[1].strip() or None, msg)])
        return _issues_out([Issue("apply", None, msg)])
    try:
        _write_atomic(ann, args.out)
    except OSError as exc:
        raise _UserError(f"cannot write {args.out}: {exc}") from None
    return 0


# -- entry point -------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ann_meta",
                                 description="Read and rewrite the metadata (COMMON / source) of a DDBJ MSS .ann")
    sub = ap.add_subparsers(dest="command", required=True)

    ex = sub.add_parser("extract", help="print the metadata as JSON (and the entry table as TSV)")
    ex.add_argument("--ann", required=True)
    ex.add_argument("--fasta", help="FASTA of the .ann (sequence lengths)")
    ex.add_argument("--json", action="store_true", help="accepted for symmetry; output is always JSON")
    ex.add_argument("--entries-out", metavar="TSV", help="also write the entry table here")

    apl = sub.add_parser("apply", help="rewrite COMMON / source from a common JSON and an entry table")
    apl.add_argument("--ann", required=True)
    apl.add_argument("--common", required=True, help="common JSON (DBLINK, SUBMITTER, ..., SOURCE, _extra)")
    apl.add_argument("--entries", required=True, help="entry table TSV (as written by extract --entries-out)")
    apl.add_argument("--out", required=True, help="output .ann (may equal --ann)")
    apl.add_argument("--category", required=True, help="submission category (GNM, WGS, ..., NONE, GENE)")
    apl.add_argument("--fasta", help="FASTA: checks seq_id and gives sequence lengths")
    apl.add_argument("--allow-unknown", metavar="Q1,Q2",
                     help="qualifier names allowed besides the known ones "
                          "(default: the unknown names already in --ann)")
    args = ap.parse_args(argv)

    try:
        return _cmd_extract(args) if args.command == "extract" else _cmd_apply(args)
    except _UserError as exc:
        _err(str(exc))
        return 1


if __name__ == "__main__":
    sys.exit(main())
