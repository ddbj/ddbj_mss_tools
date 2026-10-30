"""Data model for the metadata of one .ann file (what `ann_meta extract` reads).

``MetaDoc.to_json()`` gives plain JSON types (tuples become 2-element lists);
``MetaDoc.from_json()`` restores the dataclasses, so the two round-trip.
"""
from __future__ import annotations

import copy
from dataclasses import asdict, dataclass, field

from ann_meta.category import CategoryGuess

QualValue = str | bool  # flag qualifiers (no value) are True, as in mss_builder's JSON


@dataclass
class Entry:
    seq_id: str
    length: int | None
    type: str = "unplaced"          # chromosome | organelle | plasmid | segment | unplaced
    seq_name: str = ""
    status: str = "partial"         # complete | partial
    topology: str = ""              # circular | linear | ""
    overrides: dict[str, QualValue] = field(default_factory=dict)
    ff_definition: str | None = None
    extra: list[tuple[str, QualValue]] = field(default_factory=list)


@dataclass
class Note:
    code: str
    entry: str | None
    detail: str


@dataclass
class MetaDoc:
    category: CategoryGuess
    common: dict
    common_extra: dict[str, list[tuple[str, str]]]
    untouched_common_features: list[str]
    source: dict
    source_identifier: str | None
    entries: list[Entry]
    layout: str                     # common_source | per_entry
    notes: list[Note]
    unknown_qualifiers: list[str] = field(default_factory=list)

    def to_json(self) -> dict:
        return {
            "category": asdict(self.category),
            "common": copy.deepcopy(self.common),
            "common_extra": {k: [list(t) for t in v] for k, v in self.common_extra.items()},
            "untouched_common_features": list(self.untouched_common_features),
            "source": dict(self.source),
            "source_identifier": self.source_identifier,
            "entries": [
                {**asdict(e), "extra": [list(t) for t in e.extra]} for e in self.entries
            ],
            "layout": self.layout,
            "notes": [asdict(n) for n in self.notes],
            "unknown_qualifiers": list(self.unknown_qualifiers),
        }

    @classmethod
    def from_json(cls, d: dict) -> "MetaDoc":
        d = copy.deepcopy(d)
        cat = d.get("category") or {}
        return cls(
            category=CategoryGuess(value=cat.get("value"), evidence=list(cat.get("evidence", [])),
                                   tpa=bool(cat.get("tpa", False))),
            common=d.get("common", {}),
            common_extra={k: [(q, v) for q, v in rows] for k, rows in d.get("common_extra", {}).items()},
            untouched_common_features=list(d.get("untouched_common_features", [])),
            source=dict(d.get("source", {})),
            source_identifier=d.get("source_identifier"),
            entries=[
                Entry(**{**e, "extra": [(q, v) for q, v in e.get("extra", [])],
                         "overrides": dict(e.get("overrides", {}))})
                for e in d.get("entries", [])
            ],
            layout=d.get("layout", "per_entry"),
            notes=[Note(**n) for n in d.get("notes", [])],
            unknown_qualifiers=list(d.get("unknown_qualifiers", [])),
        )
