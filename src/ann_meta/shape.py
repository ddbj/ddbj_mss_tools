"""Structural check of the common JSON for `ann_meta apply`.

``common.models.CommonModel`` also enforces content rules (DBLINK needs a
project and a sample, SUBMITTER needs an author, the REFERENCE status enum,
``year`` as an int ...). Those are *required items*, which the metadata editor
shows but never uses as a gate, and a hand-made .ann often lacks some of them.
:class:`LenientCommon` checks only the structure: the same sections, every
field optional, values that are strings, lists of strings or (``SOURCE``)
booleans, unknown fields allowed with the same value types, no validators.
``null`` means "not given" and is dropped by :func:`normalize_common`.

:func:`normalize_common` also maps the DBLINK aliases mss_builder accepts
(``bioproject``, ``biosample``, ``sra``, ``DRA``, ``sequence read archive``)
to the MSS names before the rows are written; :func:`alias_clashes` reports a
section that gives one item under two names.
"""
from __future__ import annotations

import copy
from typing import Optional

from pydantic import BaseModel, ConfigDict, StrictBool, StrictStr

# alias -> name used by create_common (which writes sample / sequence_read_archive
# as "biosample" / "sequence read archive")
DBLINK_ALIASES = {
    "bioproject": "project",
    "biosample": "sample",
    "sequence read archive": "sequence_read_archive",
    "sra": "sequence_read_archive",
    "DRA": "sequence_read_archive",
}
# the sections apply writes (extract.COMMON_TARGETS)
_TARGETS = ("DBLINK", "SUBMITTER", "REFERENCE", "COMMENT", "DATE", "ST_COMMENT")

Str = Optional[StrictStr]
StrList = Optional[list[StrictStr]]
StrOrList = Optional[StrictStr | list[StrictStr]]
StrDict = dict[str, StrOrList]


class _Section(BaseModel):
    model_config = ConfigDict(extra="allow")
    __pydantic_extra__: dict[str, StrOrList]


class LenientDblink(_Section):
    project: Str = None
    sample: StrOrList = None
    sequence_read_archive: StrOrList = None


class LenientSubmitter(_Section):
    ab_name: StrList = None
    consrtm: Str = None
    contact: Str = None
    email: Str = None
    url: Str = None
    institute: Str = None
    department: Str = None
    country: Str = None
    state: Str = None
    city: Str = None
    street: Str = None
    zip: Str = None


class LenientReference(_Section):
    title: Str = None
    ab_name: StrList = None
    consrtm: Str = None
    status: Str = None
    year: Str = None
    journal: Str = None
    volume: Str = None
    start_page: Str = None
    end_page: Str = None


class LenientCommon(BaseModel):
    model_config = ConfigDict(extra="allow")

    DBLINK: Optional[LenientDblink] = None
    SUBMITTER: Optional[LenientSubmitter] = None
    REFERENCE: Optional[list[LenientReference]] = None
    COMMENT: Optional[StrDict | list[StrDict]] = None
    DATE: Optional[dict[str, Str]] = None
    ST_COMMENT: Optional[StrDict | list[StrDict]] = None
    SOURCE: Optional[dict[str, Optional[StrictStr | StrictBool | list[StrictStr]]]] = None
    DATATYPE: Optional[StrDict] = None
    KEYWORD: Optional[StrDict] = None
    SOURCE_IDENTIFIER: Str = None
    INFRASPECIFIC_NAME_MODIFIER: Str = None


def alias_clashes(common: dict) -> list[str]:
    """DBLINK names given twice (directly and by an alias, or by two aliases)."""
    dblink = common.get("DBLINK")
    if not isinstance(dblink, dict):
        return []
    seen: dict[str, str] = {}
    out: list[str] = []
    for key in dblink:
        name = DBLINK_ALIASES.get(key, key)
        if name in seen:
            out.append(f"DBLINK.{key}")
        else:
            seen[name] = key
    return out


def _drop_none(value):
    if isinstance(value, dict):
        return {k: _drop_none(v) for k, v in value.items() if v is not None}
    if isinstance(value, list):
        return [_drop_none(v) for v in value if v is not None]
    return value


def normalize_common(common: dict) -> dict:
    """A copy of *common* with DBLINK aliases renamed and ``null`` fields of the
    COMMON targets removed. The first of two clashing names wins (the CLI
    rejects clashes before this)."""
    out = copy.deepcopy(common)
    for name in _TARGETS:
        if out.get(name) is not None:
            out[name] = _drop_none(out[name])
    dblink = out.get("DBLINK")
    if isinstance(dblink, dict):
        renamed: dict = {}
        for key, value in dblink.items():
            renamed.setdefault(DBLINK_ALIASES.get(key, key), value)
        out["DBLINK"] = renamed
    return out
