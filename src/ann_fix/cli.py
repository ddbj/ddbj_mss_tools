"""ann_fix: list ddbj-validator findings per code and apply chosen fixes on the .ann."""
from __future__ import annotations

import argparse
import json
import sys

from .annfile import AnnFile
from .operations import ARTIFICIAL_LOCATION_VALUES, DEFAULTS, IGNORABLE, KEEP_VALUES, OPERATIONS, apply
from .report import Finding, group_by_code, parse_details


def _default_label(code: str) -> str:
    if code in DEFAULTS:
        return DEFAULTS[code][0]
    if code in IGNORABLE:
        return f"ignore ({IGNORABLE[code]})"
    return "none"


def _examples(ann: AnnFile, findings: list[Finding], n: int = 5) -> list[dict]:
    out = []
    for f in findings[:n]:
        feat = ann.feature_at_line(f.line) if f.line is not None else None
        out.append({"line": f.line, "key": feat.key if feat else f.feature,
                    "location": feat.location if feat else None,
                    "locus_tag": feat.get(ann.rows, "locus_tag") if feat else None})
    return out


def _parse_apply(spec: str, present: list[str]) -> dict[str, tuple[str, dict]]:
    if spec == "defaults":
        return {c: DEFAULTS[c] for c in present if c in DEFAULTS}
    choice: dict[str, tuple[str, dict]] = {}
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        code, eq, rest = item.partition("=")
        if not eq or not rest:
            raise ValueError(f"malformed --apply item {item!r} (expected CODE=OP[:ARG])")
        op, _, arg = rest.partition(":")
        if op not in OPERATIONS:
            raise ValueError(f"unknown operation {op!r}; choose from {sorted(OPERATIONS)}")
        kwargs = dict(DEFAULTS.get(code, ("", {}))[1]) if DEFAULTS.get(code, ("", {}))[0] == op else {}
        if arg:
            if op.startswith("overlap"):
                keep = {"trna": "other"}.get(arg.lower(), arg.lower())
                if keep not in KEEP_VALUES:
                    raise ValueError(f"keep must be one of cds/other (got {arg!r})")
                kwargs["keep"] = keep
            elif op == "artificial-location":
                if arg not in ARTIFICIAL_LOCATION_VALUES:
                    raise ValueError(f"artificial_location value must be one of {sorted(ARTIFICIAL_LOCATION_VALUES)}")
                kwargs["value"] = arg
            else:
                raise ValueError(f"{op} takes no :ARG")
        choice[code] = (op, kwargs)
    return choice


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ann_fix",
                                 description="Resolve ddbj-validator findings directly on a gff2mss .ann")
    ap.add_argument("--ann", required=True)
    ap.add_argument("--report", required=True, help="validation_report_details.txt from ddbj-validator")
    ap.add_argument("--detect", action="store_true", help="list findings per code (default when --apply is absent)")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--apply", metavar="SPEC",
                    help="'defaults' or CODE=OP[:ARG],... (OP: misc | artificial-location | drop | overlap-drop | overlap-misc)")
    ap.add_argument("--out", help="output .ann (required with --apply; may equal --ann)")
    args = ap.parse_args(argv)

    ann = AnnFile.read(args.ann)
    with open(args.report, encoding="utf-8") as fh:
        groups = group_by_code(parse_details(fh.read()))

    if args.apply is None:
        rows = [{"code": c, "level": fs[0].level, "count": len(fs), "default": _default_label(c),
                 "examples": _examples(ann, fs)} for c, fs in groups.items()]
        if args.json:
            print(json.dumps(rows, ensure_ascii=False, indent=2))
        else:
            for r in rows:
                print(f"{r['code']}  {r['level']}  n={r['count']}  default={r['default']}")
                for e in r["examples"]:
                    if e["line"] is not None:
                        print(f"  Line[{e['line']}] {e['key']} {e['location']} {e['locus_tag'] or ''}")
        return 0

    if not args.out:
        sys.stderr.write("error: --out is required with --apply\n")
        return 2
    try:
        choice = _parse_apply(args.apply, list(groups))
    except ValueError as exc:
        sys.stderr.write(f"error: {exc}\n")
        return 2
    findings = [f for fs in groups.values() for f in fs]
    try:
        log = apply(ann, findings, choice)
    except ValueError as exc:
        sys.stderr.write(f"error: {exc}; nothing written\n")
        return 1
    for line in log:
        print(line)
    ann.write(args.out)
    print(f"wrote {args.out}")
    print("re-validate: ddbj-validator ddbj <dir containing the .ann and .fasta> -o <dir>_out -f -j 1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
