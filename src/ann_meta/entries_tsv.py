"""Entry table TSV: one row per entry with its role and per-entry qualifiers.

Header: ``#seq_id type seq_name status topology ff_definition`` followed by one
column per qualifier name (the sorted union of every entry's ``overrides`` and
``extra`` names). An empty cell means "not specified". Flag qualifiers
(``True``) are written as ``yes``. A qualifier that appears more than once in
an entry's ``extra`` cannot be shown in one cell: its cell is left empty and
the full list stays in the JSON (``MetaDoc``).

``read_entries`` reads the table back: known qualifier names go to
``overrides``, unknown names and role names (``chromosome`` / ``plasmid`` /
``organelle`` / ``segment``, which ``extract`` leaves in ``extra`` when an entry
has a second role) go to ``extra`` in column order. ``yes`` becomes ``True``;
an empty cell means "not specified" (a fixed column then keeps the ``Entry``
default).
"""
from __future__ import annotations

from pathlib import Path

from ann_meta.model import Entry, QualValue
from common.qualifiers import known_qualifier_names

FIXED_COLUMNS = ("#seq_id", "type", "seq_name", "status", "topology", "ff_definition")
_ROLE_NAMES = ("chromosome", "plasmid", "organelle", "segment")


def _cell(value) -> str:
    if value is True:
        return "yes"
    if value is False or value is None:
        return ""
    return str(value)


def write_entries(entries: list[Entry], path: str | Path) -> None:
    names: set[str] = set()
    for e in entries:
        names.update(e.overrides)
        names.update(q for q, _ in e.extra)
    qual_cols = sorted(names)
    lines = ["\t".join(FIXED_COLUMNS + tuple(qual_cols))]
    for e in entries:
        counts: dict[str, int] = {}
        for q, _ in e.extra:
            counts[q] = counts.get(q, 0) + 1
        extra_single = {q: v for q, v in e.extra if counts[q] == 1}
        cells = [e.seq_id, e.type, e.seq_name, e.status, e.topology, e.ff_definition or ""]
        for q in qual_cols:
            if q in e.overrides:
                cells.append(_cell(e.overrides[q]))
            elif counts.get(q, 0) == 1:
                cells.append(_cell(extra_single[q]))
            else:
                cells.append("")
        lines.append("\t".join(cells))
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\n".join(lines) + "\n")


def _value(cell: str) -> QualValue:
    return True if cell == "yes" else cell


def read_entries(path: str | Path) -> list[Entry]:
    """Read an entry table written by :func:`write_entries` (or edited by hand).

    Raises ``ValueError`` (with the line number) for a missing / malformed header,
    duplicate column names or a row with more cells than the header.
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
    dupes = sorted({c for c in header if header.count(c) > 1})
    if dupes:
        raise ValueError(f"entry table: duplicate column(s): {', '.join(dupes)}")
    known = known_qualifier_names()
    entries: list[Entry] = []
    for lineno, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        cells = line.split("\t")
        if len(cells) > len(header):
            raise ValueError(f"entry table line {lineno}: {len(cells)} cells, header has {len(header)}")
        row = dict(zip(header, cells + [""] * (len(header) - len(cells))))
        e = Entry(seq_id=row["#seq_id"], length=None)
        for col, attr in (("type", "type"), ("seq_name", "seq_name"), ("status", "status"),
                          ("topology", "topology"), ("ff_definition", "ff_definition")):
            if row[col]:
                setattr(e, attr, row[col])
        for name in header:
            if name in FIXED_COLUMNS or not row[name]:
                continue
            if name in known and name not in _ROLE_NAMES:
                e.overrides[name] = _value(row[name])
            else:
                e.extra.append((name, _value(row[name])))
        entries.append(e)
    return entries
