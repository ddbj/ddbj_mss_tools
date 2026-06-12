"""
Source feature builders shared between egapx2mss and wgs_maker.
"""

from __future__ import annotations

from typing import Optional

Row = list[str]


# ── SequenceRoleEntry ─────────────────────────────────────────────────────────

class SequenceRoleEntry:
    """Represents one row in a sequence role file (sequence_roles.tsv; legacy name: chromosomes.txt)."""
    __slots__ = ("seq_id", "type", "seq_name", "status", "is_circular")

    def __init__(self, seq_id: str, type_: str, seq_name: str, status: str, is_circular: bool):
        self.seq_id = seq_id
        self.type = type_            # "chromosome" | "organelle" | "unplaced"
        self.seq_name = seq_name    # e.g. "1", "2", "mitochondrion", ""
        self.status = status        # "complete" | "partial"
        self.is_circular = is_circular


def load_sequence_roles(path: str) -> dict[str, SequenceRoleEntry]:
    """
    Parse a 5-column TSV sequence role file (sequence_roles.tsv; legacy name: chromosomes.txt).

    Columns: seq_id <TAB> type <TAB> seq_name <TAB> status <TAB> topology
    Lines starting with '#' are treated as header/comments and skipped.

    Returns a dict: seq_id → SequenceRoleEntry
    """
    result: dict[str, SequenceRoleEntry] = {}
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                raise ValueError(
                    f"{path}:{lineno}: expected at least 2 tab-separated columns, got {len(cols)}"
                )
            seq_id   = cols[0].strip()
            type_    = cols[1].strip() if len(cols) > 1 else "unplaced"
            seq_name = cols[2].strip() if len(cols) > 2 else ""
            status   = cols[3].strip() if len(cols) > 3 else "partial"
            topology = cols[4].strip() if len(cols) > 4 else "linear"
            is_circular = topology.lower() == "circular"
            result[seq_id] = SequenceRoleEntry(seq_id, type_, seq_name, status, is_circular)
    return result


# ── egapx2mss source qualifier helpers ───────────────────────────────────────

def source_qualifier(entry: Optional[SequenceRoleEntry], seq_id: str,
                     is_wgs: bool = False) -> dict[str, str]:
    """
    Return extra source qualifiers (beyond common.SOURCE) for one sequence entry.

    Rules follow mss_format.md:
    - unplaced + WGS: submitter_seqid = seq_id
    - unplaced + non-WGS: no extra qualifier
    - chromosome: chromosome = seq_name (omitted when seq_name is empty)
    - organelle:  organelle  = seq_name
    """
    if entry is None or entry.type == "unplaced":
        return {"submitter_seqid": seq_id} if is_wgs else {}
    if entry.type == "chromosome":
        return {"chromosome": entry.seq_name} if entry.seq_name else {}
    if entry.type == "organelle":
        return {"organelle": entry.seq_name}
    return {}


def _molecule_token(mol_type: str | None) -> str:
    """Decide the molecule token used in ff_definition ("<prefix> <token>, ...").

    Rules (in order):
      1. empty / None              -> "DNA"
      2. contains tRNA/rRNA/mRNA   -> that token (case-SENSITIVE; INSDC fixed spelling)
      3. lowercased contains "dna" -> "DNA"
      4. lowercased contains "rna" -> "RNA"
      5. otherwise (e.g. "protein")-> "DNA" (default)
    """
    if not mol_type:
        return "DNA"
    for token in ("tRNA", "rRNA", "mRNA"):
        if token in mol_type:
            return token
    low = mol_type.lower()
    if "dna" in low:
        return "DNA"
    if "rna" in low:
        return "RNA"
    return "DNA"


def ff_definition(entry: Optional[SequenceRoleEntry], seq_id: str, organism: str,
                  infraspecific_name_modifier: str, is_wgs: bool = False) -> str:
    """
    Build the ff_definition qualifier value following mss_format.md.

    *infraspecific_name_modifier* is the value of the qualifier named by SOURCE_IDENTIFIER
    (e.g. the value of 'strain' or 'isolate') from common.SOURCE.

    *is_wgs* is True when all entries in the submission are unplaced (WGS mode).
    """
    prefix = f"{organism} {infraspecific_name_modifier}".strip() if infraspecific_name_modifier else organism

    if entry is None or entry.type == "unplaced":
        if is_wgs:
            return f"{prefix} DNA, {seq_id}"
        else:
            return f"{prefix} DNA, unplaced sequence {seq_id}"

    if entry.type == "chromosome":
        chr_part = f"chromosome {entry.seq_name}".strip() if entry.seq_name else "chromosome"
        if entry.status == "complete":
            return f"{prefix} DNA, {chr_part}, complete sequence"
        else:
            return f"{prefix} DNA, {chr_part}, unlocalized sequence {seq_id}"

    if entry.type == "organelle":
        organelle_name = entry.seq_name
        if entry.status == "complete":
            return f"{prefix} DNA, {organelle_name}, complete sequence"
        else:
            return f"{prefix} DNA, {organelle_name}, partial sequence"

    # fallback
    return f"{prefix} DNA, {seq_id}"


# ── wgs_maker source feature builder ─────────────────────────────────────────

# INSDC source qualifiers written without a value (the bare key).
_FLAG_QUALIFIERS = frozenset({
    "environmental_sample",
    "transgenic",
    "germline",
    "rearranged",
    "proviral",
    "macronuclear",
    "metagenomic",
    "focus",
})

# Values interpreted as "off" for a flag qualifier (case-insensitive).
_FALSE_VALUES = frozenset({"false", "no"})


def _flag_is_set(value) -> bool:
    """Interpret a flag-qualifier value as on/off.

    Off: boolean False, or the strings "false"/"no" (case-insensitive,
    surrounding whitespace ignored).
    On:  everything else — boolean True, "" (empty string, kept on for
    backward compatibility with ``"environmental_sample": ""``), "true",
    "yes", "1", and any other non-empty string.
    """
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() not in _FALSE_VALUES


def _source_qualifier_rows(key: str, value) -> list[Row]:
    """Emit source qualifier row(s) for one (key, value).

    Flag qualifiers (in ``_FLAG_QUALIFIERS``) emit a single valueless row
    when the value is truthy (see :func:`_flag_is_set`) and are omitted
    entirely when falsy. Other qualifiers emit one row per value — a list
    yields one row per element (e.g. multiple ``culture_collection``).
    """
    if key in _FLAG_QUALIFIERS:
        return [["", "", "", key, ""]] if _flag_is_set(value) else []
    if isinstance(value, list):
        return [["", "", "", key, str(v)] for v in value]
    return [["", "", "", key, str(value)]]


def create_source_feature(
    _submission_category: str,
    seq_name: Optional[str],
    seq_type: Optional[str],
    seq_topology: Optional[str],
    source_dict: dict,
    source_modifier_key: str = "",
    use_meta_expression: bool = False,
) -> list[Row]:
    """
    Build source feature rows for DDBJ MSS annotation.

    Returns a list of 5-column rows.

    When *use_meta_expression* is True the feature is built for inclusion in the
    COMMON block: location is ``1..E``, ``submitter_seqid`` is ``@@[entry]@@``,
    and ``ff_definition`` uses ``@@[...]@@`` meta-notation.  *source_modifier_key*
    (from ``SOURCE_IDENTIFIER``, e.g. ``"strain"``, ``"cultivar"``) is used
    to select the modifier placeholder in ``ff_definition``; when empty only
    ``@@[organism]@@`` is included.
    """
    if use_meta_expression:
        return _create_source_with_meta(
            _submission_category, source_dict, source_modifier_key
        )

    from common.submission_category import get_category_rules
    rules = get_category_rules(_submission_category)
    env_sample = "environmental_sample" in rules.auto_source_qualifiers
    modifier = "isolate" if env_sample else "strain"

    submitter_seqid = None
    plasmid = False
    mol_type = source_dict.get("mol_type", "genomic DNA")

    # WGS-family: source goes in COMMON block, submitter_seqid always set
    if rules.datatype == "WGS":
        submitter_seqid = "@@[entry]@@"
        ff_def = f"@@[organism]@@ @@[{modifier}]@@ DNA, @@[submitter_seqid]@@"
    else:
        # Per-entry source (GNM, MAG, etc.)
        if seq_type in ["c", "complete"]:
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ DNA, complete genome"
        elif seq_type in ["n", "nearly complete", "nearly-complete"]:
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ DNA, nearly complete genome"
        elif seq_type in ["p", "plasmid"]:
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ plasmid @@[plasmid]@@ DNA, complete sequence"
            plasmid = True
        else:
            submitter_seqid = "@@[entry]@@"
            ff_def = f"@@[organism]@@ @@[{modifier}]@@ DNA, @@[submitter_seqid]@@"

    ret: list[Row] = []
    ret.append(["", "source", "1..E", "mol_type", mol_type])
    ret.append(["", "", "", "ff_definition", ff_def])
    if submitter_seqid:
        ret.append(["", "", "", "submitter_seqid", submitter_seqid])
    if env_sample:
        ret.append(["", "", "", "environmental_sample", ""])
    if plasmid:
        ret.append(["", "", "", "plasmid", seq_name])
    auto_keys = set(rules.auto_source_qualifiers.keys())
    for key, value in source_dict.items():
        if key != "mol_type" and key not in auto_keys:
            ret.extend(_source_qualifier_rows(key, value))
    if rules.datatype == "WGS":
        # Source feature will be appended to COMMON. Nothing more to do.
        pass
    else:
        if seq_topology in ["c", "circular"]:
            ret = [["", "TOPOLOGY", "", "circular", ""]] + ret
        ret[0][0] = seq_name
    return ret


def _create_source_with_meta(
    category: str,
    source_dict: dict,
    source_modifier_key: str,
) -> list[Row]:
    """Build source feature rows using ``@@[...]@@`` meta-notation for the COMMON block.

    Location is ``1..E`` (expands to the full length of each entry).
    ``submitter_seqid`` is set to ``@@[entry]@@`` so DDBJ MSS substitutes the
    per-entry identifier at submission time.
    ``ff_definition`` references ``@@[organism]@@`` and, when *source_modifier_key*
    is provided, ``@@[<source_modifier_key>]@@``.
    """
    from common.submission_category import get_category_rules
    rules = get_category_rules(category)
    environmental_sample = "environmental_sample" in rules.auto_source_qualifiers
    mol_type = source_dict.get("mol_type", "genomic DNA")

    if source_modifier_key:
        ff_def = f"@@[organism]@@ @@[{source_modifier_key}]@@ DNA, @@[submitter_seqid]@@"
    else:
        ff_def = "@@[organism]@@ DNA, @@[submitter_seqid]@@"

    rows: list[Row] = []
    rows.append(["", "source", "1..E", "mol_type", mol_type])
    rows.append(["", "", "", "ff_definition", ff_def])
    rows.append(["", "", "", "submitter_seqid", "@@[entry]@@"])
    if environmental_sample:
        rows.append(["", "", "", "environmental_sample", ""])
    auto_keys = set(rules.auto_source_qualifiers.keys())
    for key, value in source_dict.items():
        if key != "mol_type" and key not in auto_keys:
            rows.extend(_source_qualifier_rows(key, value))
    return rows
