"""Guess the DDBJ MSS submission data type ("category") of an .ann file.

The guess is made from plain data -- the COMMON block, each entry's source
qualifiers, whether an entry has a TOPOLOGY row, its ff_definition, and the
first entry's length -- so that callers can feed it from any parser.

Rules are the table RULES below, tried from top to bottom; the first one
with any evidence wins (TPA is a separate flag and is checked before
everything else). Source: Confluence 2026-07-16 MSS Meeting.

    TLS  DATATYPE TLS, or KEYWORD TLS / Targeted Locus Study
    TSA  DATATYPE TSA, DIVISION TSA, or KEYWORD TSA / Transcriptome Shotgun Assembly
    HTG  DATATYPE HTG, DIVISION HTG, or KEYWORD HTG / HTGS_PHASE* / HTGS_*
    HTC  DIVISION HTC, or KEYWORD HTC / HTC_FLI / FLI_CDNA
    EST  DIVISION EST, or KEYWORD EST
    MAG  KEYWORD MAG / Metagenome Assembled Genome, or DIVISION ENV with a
         source /metagenome_source
    SAG  KEYWORD SAG / Single Amplified Genome
    WGS  DATATYPE WGS (or KEYWORD WGS) not caught by MAG / SAG above
    GNM  KEYWORD GNM / Finished level genomic sequences
    GNM  no DATATYPE, an organism contains "virus" / "phage", and BioProject
         and BioSample are both given (COMMON DBLINK or that entry's DBLINK)
    MISC no DATATYPE, and every ff_definition is "plasmid ... complete
         sequence", or every entry has /organelle, or an organism contains
         "virus" / "phage" (a virus/phage without BioProject+BioSample lands
         here; DDBJ files plasmid-only and organelle-only submissions under
         MISC even though they are complete and usually carry TOPOLOGY)
    GNM  no DATATYPE, and a TOPOLOGY row, an ff_definition with
         "complete genome" / ", chromosome", or a first entry > 1 Mb
    None nothing above matched (evidence is empty)

Keyword matching is case-insensitive on whole keyword values; a KEYWORD
value may itself hold several keywords separated by ";".

Whether a MAG is a draft or a complete genome (DATATYPE WGS or not) is left
to the caller.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

# mss-reporter filetype.py COMPLETE_GENOME_MIN_LENGTH
COMPLETE_GENOME_MIN_LENGTH = 1_000_000

# Allowed CategoryGuess.value strings, for callers (form choices, validation) in later ann_meta tasks.
CATEGORIES = ("GNM", "WGS", "MAG", "SAG", "TLS", "HTG", "TSA", "HTC", "EST", "MISC")


@dataclass
class EntrySketch:
    seq_id: str
    has_topology: bool
    ff_definition: str | None
    qualifiers: dict[str, str] = field(default_factory=dict)  # source qualifiers
    dblink: dict[str, list[str]] = field(default_factory=dict)  # the entry's own DBLINK, {mss_name: [values]}


@dataclass
class CategoryGuess:
    value: str | None
    evidence: list[str] = field(default_factory=list)
    tpa: bool = False


def _values(common: dict[str, list[tuple[str, str]]], feature: str, qualifier: str) -> list[str]:
    return [v.strip() for q, v in common.get(feature, []) if q == qualifier and v is not None and v.strip()]


def _keywords(common: dict[str, list[tuple[str, str]]]) -> list[str]:
    out: list[str] = []
    for _q, v in common.get("KEYWORD", []):
        if v is None:
            continue
        out.extend(k.strip() for k in v.split(";") if k.strip())
    return out


def _kw_evidence(keywords: list[str], exact: set[str], prefixes: tuple[str, ...] = ()) -> list[str]:
    return [f"KEYWORD {k}" for k in keywords
            if k.upper() in exact or any(k.upper().startswith(p) for p in prefixes)]


def _eq_evidence(values: list[str], label: str, wanted: str) -> list[str]:
    return [label.format(v) for v in values if v.upper() == wanted]


_DBLINK_ALIASES = {"project": "project", "bioproject": "project", "biosample": "sample", "sample": "sample"}


@dataclass
class _Facts:
    common: dict[str, list[tuple[str, str]]]
    entries: list[EntrySketch]
    first_entry_length: int | None
    datatypes: list[str]
    divisions: list[str]
    keywords: list[str]

    def dt(self, wanted: str) -> list[str]:
        return _eq_evidence(self.datatypes, "DATATYPE type={}", wanted)

    def div(self, wanted: str) -> list[str]:
        return _eq_evidence(self.divisions, "DIVISION division={}", wanted)

    def kw(self, exact: set[str], prefixes: tuple[str, ...] = ()) -> list[str]:
        return _kw_evidence(self.keywords, exact, prefixes)

    def defs(self) -> list[str]:
        return [(e.ff_definition or "").lower() for e in self.entries]


@dataclass(frozen=True)
class Rule:
    value: str
    test: Callable[[_Facts], list[str]]


def _tls(f: _Facts) -> list[str]:
    return f.dt("TLS") + f.kw({"TLS", "TARGETED LOCUS STUDY"})


def _tsa(f: _Facts) -> list[str]:
    return f.dt("TSA") + f.div("TSA") + f.kw({"TSA", "TRANSCRIPTOME SHOTGUN ASSEMBLY"})


def _htg(f: _Facts) -> list[str]:
    return f.dt("HTG") + f.div("HTG") + f.kw({"HTG"}, ("HTGS_",))


def _htc(f: _Facts) -> list[str]:
    return f.div("HTC") + f.kw({"HTC", "HTC_FLI", "FLI_CDNA"})


def _est(f: _Facts) -> list[str]:
    return f.div("EST") + f.kw({"EST"})


def _mag(f: _Facts) -> list[str]:
    # ddbj-validator definitions.json: MAG required_keywords [["MAG"], ["Metagenome Assembled Genome"]]
    out = f.kw({"MAG", "METAGENOME ASSEMBLED GENOME"})
    env = f.div("ENV")
    if env and any("metagenome_source" in e.qualifiers for e in f.entries):
        out += env + ["source metagenome_source"]
    return out


def _sag(f: _Facts) -> list[str]:
    # ddbj-validator definitions.json:6540 defines SAG without any keyword rule
    # (only isolation_source, ST_COMMENT and invalid_moltype), so SAG is
    # detected from KEYWORD alone.
    return f.kw({"SAG", "SINGLE AMPLIFIED GENOME"})


def _wgs(f: _Facts) -> list[str]:
    return f.dt("WGS") + f.kw({"WGS"})


def _gnm_keyword(f: _Facts) -> list[str]:
    return f.kw({"GNM", "FINISHED LEVEL GENOMIC SEQUENCES"})


def _dblink_names(common: dict[str, list[tuple[str, str]]], entry: EntrySketch) -> set[str]:
    names = {_DBLINK_ALIASES[q.strip().lower()] for q, v in common.get("DBLINK", [])
             if q.strip().lower() in _DBLINK_ALIASES and v is not None and v.strip()}
    for q, vals in entry.dblink.items():
        name = _DBLINK_ALIASES.get(q.strip().lower())
        if name and any(v and v.strip() for v in vals):
            names.add(name)
    return names


def _gnm_virus(f: _Facts) -> list[str]:
    if f.datatypes:
        return []
    for e in f.entries:
        if any(w in (e.qualifiers.get("organism") or "").lower() for w in ("virus", "phage")) \
                and {"project", "sample"} <= _dblink_names(f.common, e):
            word = "virus" if "virus" in e.qualifiers["organism"].lower() else "phage"
            return [f"organism contains '{word}'", "DBLINK project and biosample"]
    return []


def _misc(f: _Facts) -> list[str]:
    if f.datatypes:
        return []
    defs = f.defs()
    out: list[str] = []
    if f.entries and all("plasmid" in d and "complete sequence" in d for d in defs):
        out.append("ff_definition 'plasmid ... complete sequence' (all entries)")
    if f.entries and all("organelle" in e.qualifiers for e in f.entries):
        out.append("source organelle (all entries)")
    organisms = [(e.qualifiers.get("organism") or "").lower() for e in f.entries]
    for word in ("virus", "phage"):
        if any(word in o for o in organisms):
            out.append(f"organism contains '{word}'")
    return out


def _gnm_structure(f: _Facts) -> list[str]:
    if f.datatypes:
        return []
    defs = f.defs()
    out: list[str] = []
    if any(e.has_topology for e in f.entries):
        out.append("TOPOLOGY present")
    for phrase in ("complete genome", ", chromosome"):
        if any(phrase in d for d in defs):
            out.append(f"ff_definition '{phrase}'")
    n = f.first_entry_length
    if n is not None and n > COMPLETE_GENOME_MIN_LENGTH:
        out.append(f"first entry {n:,} bp > {COMPLETE_GENOME_MIN_LENGTH:,}")
    return out


RULES: tuple[Rule, ...] = (
    Rule("TLS", _tls),
    Rule("TSA", _tsa),
    Rule("HTG", _htg),
    Rule("HTC", _htc),
    Rule("EST", _est),
    Rule("MAG", _mag),
    Rule("SAG", _sag),
    Rule("WGS", _wgs),
    Rule("GNM", _gnm_keyword),
    Rule("GNM", _gnm_virus),
    Rule("MISC", _misc),
    Rule("GNM", _gnm_structure),
)


def detect_category(common: dict[str, list[tuple[str, str]]], entries: list[EntrySketch],
                    first_entry_length: int | None) -> CategoryGuess:
    facts = _Facts(common, entries, first_entry_length, _values(common, "DATATYPE", "type"),
                   _values(common, "DIVISION", "division"), _keywords(common))
    tpa = [f"DATATYPE type={v}" for v in facts.datatypes if v.upper().startswith("TPA")]
    if tpa:
        return CategoryGuess(value=None, evidence=tpa, tpa=True)
    for rule in RULES:
        evidence = rule.test(facts)
        if evidence:
            return CategoryGuess(value=rule.value, evidence=evidence, tpa=False)
    return CategoryGuess(value=None, evidence=[], tpa=False)
