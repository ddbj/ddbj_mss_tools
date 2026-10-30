"""Entry table TSV: one row per entry with its role and per-entry qualifiers.

Header: ``#seq_id type seq_name status topology ff_definition`` followed by the
qualifier columns: every name in an entry's ``overrides`` or ``extra``, sorted,
and repeated as often as the entry with the most occurrences has it (two
``/note`` rows -> two ``note`` columns; the n-th column holds the n-th
occurrence, entries with fewer leave the rest empty). The six fixed columns
are unique. An empty cell means "not specified". Flag qualifiers (``True``)
are written as ``yes``.

``read_entries`` reads the table back, row by row: a qualifier with one
non-empty cell goes to ``overrides`` when its name is known (and not a role
name), otherwise to ``extra``; a qualifier with two or more non-empty cells
(a repeated column) goes to ``extra``, all values in column order. Role names
(``chromosome`` / ``plasmid`` / ``organelle`` / ``segment``) always go to
``extra`` (``validate`` reports them: the role comes from ``type``). An empty
fixed cell keeps the ``Entry`` default.

``yes`` becomes the flag ``True`` only for a known boolean qualifier
(``environmental_sample`` ...) or an unknown name; a known non-boolean
qualifier keeps the string (``/strain="yes"``). Ambiguity: for an unknown
qualifier the table cannot tell a flag from the literal value ``yes``.
"""
from __future__ import annotations

from pathlib import Path

from ann_meta.model import Entry, QualValue
from common.qualifiers import known_qualifier_names, qualifier_def

FIXED_COLUMNS = ("#seq_id", "type", "seq_name", "status", "topology", "ff_definition")
_ROLE_NAMES = ("chromosome", "plasmid", "organelle", "segment")


def _cell(value) -> str:
    if value is True:
        return "yes"
    if value is False or value is None:
        return ""
    return str(value)


def _qualifier_values(e: Entry) -> dict[str, list[QualValue]]:
    out: dict[str, list[QualValue]] = {}
    for q, v in e.overrides.items():
        out.setdefault(q, []).append(v)
    for q, v in e.extra:
        out.setdefault(q, []).append(v)
    return out


def write_entries(entries: list[Entry], path: str | Path) -> None:
    values = [_qualifier_values(e) for e in entries]
    width: dict[str, int] = {}
    for vals in values:
        for q, vs in vals.items():
            width[q] = max(width.get(q, 0), len(vs))
    names = sorted(width)
    header = list(FIXED_COLUMNS) + [q for q in names for _ in range(width[q])]
    lines = ["\t".join(header)]
    for e, vals in zip(entries, values):
        cells = [e.seq_id, e.type, e.seq_name, e.status, e.topology, e.ff_definition or ""]
        for q in names:
            vs = vals.get(q, [])
            cells.extend(_cell(v) for v in vs)
            cells.extend([""] * (width[q] - len(vs)))
        lines.append("\t".join(cells))
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\n".join(lines) + "\n")


def _value(name: str, cell: str) -> QualValue:
    if cell != "yes":
        return cell
    d = qualifier_def(name)
    return True if d is None or d.get("boolean") else cell


def read_entries(path: str | Path) -> list[Entry]:
    """Read an entry table written by :func:`write_entries` (or edited by hand).

    Raises ``ValueError`` (with the line number) for a missing / malformed header,
    a repeated fixed column or a row with more cells than the header.
    """
    with open(path, encoding="utf-8", newline="") as fh:
        # split on LF only: str.splitlines would also split on \x1c, \u2028 ... inside a cell
        lines = [ln[:-1] if ln.endswith("\r") else ln for ln in fh.read().split("\n")]
    if not lines or not lines[0].startswith("#seq_id"):
        raise ValueError("entry table: the first line must be the header starting with '#seq_id'")
    header = lines[0].split("\t")
    missing = [c for c in FIXED_COLUMNS if c not in header]
    if missing:
        raise ValueError(f"entry table: header lacks column(s): {', '.join(missing)}")
    dupes = [c for c in FIXED_COLUMNS if header.count(c) > 1]
    if dupes:
        raise ValueError(f"entry table: duplicate column(s): {', '.join(dupes)}")
    fixed = {c: header.index(c) for c in FIXED_COLUMNS}
    qual_cols = [(i, name) for i, name in enumerate(header) if name not in FIXED_COLUMNS]
    known = known_qualifier_names()
    entries: list[Entry] = []
    for lineno, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        cells = line.split("\t")
        if len(cells) > len(header):
            raise ValueError(f"entry table line {lineno}: {len(cells)} cells, header has {len(header)}")
        cells += [""] * (len(header) - len(cells))
        e = Entry(seq_id=cells[fixed["#seq_id"]], length=None)
        for col in ("type", "seq_name", "status", "topology", "ff_definition"):
            if cells[fixed[col]]:
                setattr(e, col, cells[fixed[col]])
        given: dict[str, list[str]] = {}
        for i, name in qual_cols:
            if cells[i]:
                given.setdefault(name, []).append(cells[i])
        for name, vs in given.items():
            if len(vs) == 1 and name in known and name not in _ROLE_NAMES:
                e.overrides[name] = _value(name, vs[0])
            else:
                e.extra.extend((name, _value(name, v)) for v in vs)
        entries.append(e)
    return entries
