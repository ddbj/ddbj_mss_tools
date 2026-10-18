"""Parse ddbj-validator validation_report_details.txt lines.

    CODE:LEVEL:ENTRY:FEATURE[:QUALIFIER]:Line[N]:message
    CODE:LEVEL:ALL:ALL[:QUALIFIER]:[FILE]:message          (file-level, no line)
"""
from __future__ import annotations

import re
from dataclasses import dataclass

_LINE = re.compile(
    r"^(?P<code>[A-Z]+\d+):(?P<level>ERR|WAR|INFO?):(?P<entry>[^:]*):(?P<feature>[^:]*):"
    r"(?:(?P<qual>[^:\[]*):)?(?:Line\[(?P<line>\d+)\]|\[FILE\]|[^:]*):(?P<message>.*)$")


@dataclass(frozen=True)
class Finding:
    code: str
    level: str          # "ERR" | "WAR" | "INF"
    entry: str
    feature: str
    line: int | None    # 1-based line of the feature's first row in the .ann
    message: str


def parse_details(text: str) -> list[Finding]:
    out: list[Finding] = []
    for raw in text.splitlines():
        m = _LINE.match(raw.strip())
        if not m:
            continue
        line = m.group("line")
        out.append(Finding(m.group("code"), m.group("level")[:3], m.group("entry"), m.group("feature"),
                           int(line) if line else None, m.group("message").strip()))
    return out


def group_by_code(findings: list[Finding]) -> dict[str, list[Finding]]:
    groups: dict[str, list[Finding]] = {}
    for f in findings:
        groups.setdefault(f.code, []).append(f)
    return groups
