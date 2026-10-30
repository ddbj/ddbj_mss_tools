"""Entry table TSV: one row per entry with its role and per-entry qualifiers.

Header: ``#seq_id type seq_name status topology ff_definition`` followed by one
column per qualifier name (the sorted union of every entry's ``overrides`` and
``extra`` names). An empty cell means "not specified". Flag qualifiers
(``True``) are written as ``yes``. A qualifier that appears more than once in
an entry's ``extra`` cannot be shown in one cell: its cell is left empty and
the full list stays in the JSON (``MetaDoc``).
"""
from __future__ import annotations

from pathlib import Path

from ann_meta.model import Entry

FIXED_COLUMNS = ("#seq_id", "type", "seq_name", "status", "topology", "ff_definition")


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
