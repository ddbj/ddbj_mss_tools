"""Guess the DDBJ MSS submission data type ("category") of an .ann file.

The guess is made from plain data -- the COMMON block, each entry's source
qualifiers, whether an entry has a TOPOLOGY row, its ff_definition, and the
first entry's length -- so that callers can feed it from any parser.

Rules are tried from top to bottom and the first category with any evidence
wins (TPA is a separate flag and is checked before everything else):

    TLS  DATATYPE TLS, or KEYWORD TLS / Targeted Locus Study
    TSA  DATATYPE TSA, DIVISION TSA, or KEYWORD TSA / Transcriptome Shotgun Assembly
    HTG  DATATYPE HTG, DIVISION HTG, or KEYWORD HTG / HTGS_PHASE* / HTGS_*
    HTC  DIVISION HTC, or KEYWORD HTC / HTC_FLI / FLI_CDNA
    EST  DIVISION EST, or KEYWORD EST
    MAG  KEYWORD MAG / Metagenome Assembled Genome, or DIVISION ENV with a
         source /metagenome_source
    SAG  KEYWORD SAG / Single Amplified Genome
    WGS  DATATYPE WGS (or KEYWORD WGS) not caught by MAG / SAG above
    MISC no DATATYPE, and every ff_definition is "plasmid ... complete
         sequence", or every entry has /organelle, or an organism contains
         "virus" / "phage" (checked before GNM: DDBJ files virus/phage
         genomes, plasmid-only and organelle-only submissions under MISC
         even though they are complete and usually carry TOPOLOGY)
    GNM  no DATATYPE, and a TOPOLOGY row, an ff_definition with
         "complete genome" / "chromosome", or a first entry > 1 Mb
    None nothing above matched (evidence is empty)

Keyword matching is case-insensitive on whole keyword values; a KEYWORD
value may itself hold several keywords separated by ";".

Whether a MAG is a draft or a complete genome (DATATYPE WGS or not) is left
to the caller.
"""
from __future__ import annotations

from dataclasses import dataclass, field

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


def detect_category(common: dict[str, list[tuple[str, str]]], entries: list[EntrySketch],
                    first_entry_length: int | None) -> CategoryGuess:
    datatypes = _values(common, "DATATYPE", "type")
    divisions = _values(common, "DIVISION", "division")
    keywords = _keywords(common)

    def dt(wanted: str) -> list[str]:
        return _eq_evidence(datatypes, "DATATYPE type={}", wanted)

    def div(wanted: str) -> list[str]:
        return _eq_evidence(divisions, "DIVISION division={}", wanted)

    tpa = [f"DATATYPE type={v}" for v in datatypes if v.upper().startswith("TPA")]
    if tpa:
        return CategoryGuess(value=None, evidence=tpa, tpa=True)

    rules: list[tuple[str, list[str]]] = []
    rules.append(("TLS", dt("TLS") + _kw_evidence(keywords, {"TLS", "TARGETED LOCUS STUDY"})))
    rules.append(("TSA", dt("TSA") + div("TSA")
                  + _kw_evidence(keywords, {"TSA", "TRANSCRIPTOME SHOTGUN ASSEMBLY"})))
    rules.append(("HTG", dt("HTG") + div("HTG") + _kw_evidence(keywords, {"HTG"}, ("HTGS_",))))
    rules.append(("HTC", div("HTC") + _kw_evidence(keywords, {"HTC", "HTC_FLI", "FLI_CDNA"})))
    rules.append(("EST", div("EST") + _kw_evidence(keywords, {"EST"})))

    # ddbj-validator definitions.json: MAG required_keywords [["MAG"], ["Metagenome Assembled Genome"]]
    mag = _kw_evidence(keywords, {"MAG", "METAGENOME ASSEMBLED GENOME"})
    env = div("ENV")
    if env and any("metagenome_source" in e.qualifiers for e in entries):
        mag += env + ["source metagenome_source"]
    rules.append(("MAG", mag))

    # ddbj-validator definitions.json:6540 defines SAG without any keyword rule
    # (only isolation_source, ST_COMMENT and invalid_moltype), so SAG is
    # detected from KEYWORD alone.
    rules.append(("SAG", _kw_evidence(keywords, {"SAG", "SINGLE AMPLIFIED GENOME"})))
    rules.append(("WGS", dt("WGS") + _kw_evidence(keywords, {"WGS"})))

    if not datatypes:
        defs = [(e.ff_definition or "").lower() for e in entries]
        misc: list[str] = []
        if entries and all("plasmid" in d and "complete sequence" in d for d in defs):
            misc.append("ff_definition 'plasmid ... complete sequence' (all entries)")
        if entries and all("organelle" in e.qualifiers for e in entries):
            misc.append("source organelle (all entries)")
        organisms = [(e.qualifiers.get("organism") or "").lower() for e in entries]
        for word in ("virus", "phage"):
            if any(word in o for o in organisms):
                misc.append(f"organism contains '{word}'")
        rules.append(("MISC", misc))

        gnm: list[str] = []
        if any(e.has_topology for e in entries):
            gnm.append("TOPOLOGY present")
        for phrase in ("complete genome", "chromosome"):
            if any(phrase in d for d in defs):
                gnm.append(f"ff_definition '{phrase}'")
        if first_entry_length is not None and first_entry_length > COMPLETE_GENOME_MIN_LENGTH:
            gnm.append(f"first entry {first_entry_length:,} bp > {COMPLETE_GENOME_MIN_LENGTH:,}")
        rules.append(("GNM", gnm))

    for value, evidence in rules:
        if evidence:
            return CategoryGuess(value=value, evidence=evidence, tpa=False)
    return CategoryGuess(value=None, evidence=[], tpa=False)
