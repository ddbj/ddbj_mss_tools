"""Line-exact model of a DDBJ MSS .ann file.

Every physical line is kept as its 5 tab-separated columns plus its own line
ending, so rows no operation touches are written back byte-for-byte. A feature
is the contiguous run of rows starting at a row with a non-empty column 2.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

_NCOLS = 5


@dataclass
class AnnRow:
    cols: list[str]
    eol: str


@dataclass
class AnnFeature:
    entry: str
    key: str
    location: str
    first: int
    last: int

    def qualifiers(self, rows: list[AnnRow]) -> list[tuple[str, str]]:
        out = []
        for r in rows[self.first:self.last + 1]:
            if r.cols[3]:
                out.append((r.cols[3], r.cols[4]))
        return out

    def get(self, rows: list[AnnRow], key: str) -> str | None:
        for k, v in self.qualifiers(rows):
            if k == key:
                return v
        return None


def _split(line: str) -> AnnRow:
    if line.endswith("\r\n"):
        body, eol = line[:-2], "\r\n"
    elif line.endswith("\n"):
        body, eol = line[:-1], "\n"
    else:
        body, eol = line, ""
    cols = body.split("\t")
    cols = (cols + [""] * _NCOLS)[:_NCOLS] if len(cols) <= _NCOLS else cols
    return AnnRow(cols, eol)


class AnnFile:
    def __init__(self, rows: list[AnnRow]):
        self.rows = rows

    @classmethod
    def read(cls, path: str | Path) -> "AnnFile":
        with open(path, newline="", encoding="utf-8") as fh:
            text = fh.read()
        return cls([_split(line) for line in text.splitlines(keepends=True)])

    def write(self, path: str | Path) -> None:
        with open(path, "w", newline="", encoding="utf-8") as fh:
            for r in self.rows:
                fh.write("\t".join(r.cols) + r.eol)

    # -- structure ---------------------------------------------------------

    def features(self) -> list[AnnFeature]:
        feats: list[AnnFeature] = []
        entry = ""
        for i, r in enumerate(self.rows):
            if r.cols[0]:
                entry = r.cols[0]
            if r.cols[1]:
                if feats:
                    feats[-1].last = i - 1
                feats.append(AnnFeature(entry, r.cols[1], r.cols[2], i, i))
        if feats:
            feats[-1].last = self._block_end(feats[-1].first)
        return feats

    def _block_end(self, first: int) -> int:
        i = first
        while i + 1 < len(self.rows) and not self.rows[i + 1].cols[1] and _is_qualifier_row(self.rows[i + 1]):
            i += 1
        return i

    def feature_at_line(self, line_1based: int) -> AnnFeature | None:
        idx = line_1based - 1
        for f in self.features():
            if f.first == idx:
                return f
        return None

    # -- edits -------------------------------------------------------------

    def _entry_id_of(self, feat: AnnFeature) -> str:
        return self.rows[feat.first].cols[0]

    def default_eol(self) -> str:
        """The first non-empty line ending in the file ("\n" if there is none).

        Used where the neighbouring row has no line ending (a final line without
        a newline), so a CRLF file never gets a bare LF.
        """
        for r in self.rows:
            if r.eol:
                return r.eol
        return "\n"

    def replace_feature(self, feat: AnnFeature, key: str, location: str,
                        qualifiers: list[tuple[str, str]]) -> None:
        entry_col = self._entry_id_of(feat)
        eol = self.rows[feat.first].eol or self.default_eol()
        quals = qualifiers or [("", "")]
        new_rows = []
        for i, (k, v) in enumerate(quals):
            new_rows.append(AnnRow([entry_col if i == 0 else "", key if i == 0 else "",
                                    location if i == 0 else "", k, v], eol))
        last_eol = self.rows[feat.last].eol
        new_rows[-1].eol = last_eol if feat.last == len(self.rows) - 1 else eol
        self.rows[feat.first:feat.last + 1] = new_rows

    def delete_feature(self, feat: AnnFeature) -> None:
        entry_col = self._entry_id_of(feat)
        was_last = feat.last == len(self.rows) - 1
        last_eol = self.rows[feat.last].eol
        del self.rows[feat.first:feat.last + 1]
        if entry_col and feat.first < len(self.rows) and not self.rows[feat.first].cols[0]:
            self.rows[feat.first].cols[0] = entry_col
        if was_last and self.rows:
            self.rows[-1].eol = last_eol

    def add_qualifier(self, feat: AnnFeature, key: str, value: str) -> None:
        eol = self.rows[feat.last].eol          # "" only on a final line without newline
        if not eol:
            self.rows[feat.last].eol = self.default_eol()
        self.rows.insert(feat.last + 1, AnnRow(["", "", "", key, value], eol))

    def insert_feature(self, at_row: int, entry_col: str, key: str, location: str,
                       qualifiers: list[tuple[str, str]]) -> AnnFeature:
        """Insert a new feature before row index *at_row* (``len(rows)`` appends).

        New rows take the line ending of the row before *at_row* (at the top, of the
        row at *at_row*; when that is empty, the file's first line ending, else "\n").
        If the row at *at_row* carries *entry_col* as its entry marker, the marker
        moves to the inserted block so an entry never has two markers (another
        entry's marker is left alone).
        Appending after a final row without a newline gives that row the file's
        line ending and leaves the new last row without one, so the file keeps
        its shape.
        """
        if not 0 <= at_row <= len(self.rows):
            raise IndexError(f"insert position out of range: {at_row}")
        prev = self.rows[at_row - 1] if at_row > 0 else None
        nxt = self.rows[at_row] if at_row < len(self.rows) else None
        eol = (prev.eol if prev else (nxt.eol if nxt else "")) or self.default_eol()
        if entry_col and nxt is not None and nxt.cols[0] == entry_col:
            nxt.cols[0] = ""                    # the inserted block takes over the entry marker
        quals = qualifiers or [("", "")]
        new_rows = [AnnRow([entry_col if i == 0 else "", key if i == 0 else "",
                            location if i == 0 else "", k, v], eol)
                    for i, (k, v) in enumerate(quals)]
        if at_row == len(self.rows) and prev is not None and prev.eol == "":
            prev.eol = eol
            new_rows[-1].eol = ""
        self.rows[at_row:at_row] = new_rows
        entry = entry_col
        if not entry:
            for r in reversed(self.rows[:at_row]):
                if r.cols[0]:
                    entry = r.cols[0]
                    break
        return AnnFeature(entry, key, location, at_row, at_row + len(new_rows) - 1)


def _is_qualifier_row(r: AnnRow) -> bool:
    return not r.cols[0] and not r.cols[1] and not r.cols[2]

