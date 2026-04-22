"""
Source feature builders shared between egapx2mss and wgs_maker.
"""

from __future__ import annotations

from typing import Optional

Row = list[str]


# ── ChromosomeEntry (egapx2mss) ───────────────────────────────────────────────

class ChromosomeEntry:
    """Represents one row in chromosomes.txt."""
    __slots__ = ("seq_id", "type", "seq_name", "status", "is_circular")

    def __init__(self, seq_id: str, type_: str, seq_name: str, status: str, is_circular: bool):
        self.seq_id = seq_id
        self.type = type_            # "chromosome" | "organelle" | "unplaced"
        self.seq_name = seq_name    # e.g. "1", "2", "mitochondrion", ""
        self.status = status        # "complete" | "partial"
        self.is_circular = is_circular


def load_chromosomes(path: str) -> dict[str, ChromosomeEntry]:
    """
    Parse a 5-column TSV chromosomes.txt file.

    Columns: seq_id <TAB> type <TAB> seq_name <TAB> status <TAB> topology
    Lines starting with '#' are treated as header/comments and skipped.

    Returns a dict: seq_id → ChromosomeEntry
    """
    result: dict[str, ChromosomeEntry] = {}
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
            result[seq_id] = ChromosomeEntry(seq_id, type_, seq_name, status, is_circular)
    return result


# ── egapx2mss source qualifier helpers ───────────────────────────────────────

def source_qualifier(entry: Optional[ChromosomeEntry], seq_id: str,
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


def ff_definition(entry: Optional[ChromosomeEntry], seq_id: str, organism: str,
                  infraspecific_name_modifier: str, is_wgs: bool = False) -> str:
    """
    Build the ff_definition qualifier value following mss_format.md.

    *infraspecific_name_modifier* is the value of the qualifier named by INFRASPECIFIC_NAME_MODIFIER
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

def create_source_feature(
    _trad_submission_category: str,
    seq_name: Optional[str],
    seq_type: Optional[str],
    seq_topology: Optional[str],
    source_dict: dict,
    source_modifier_key: str = "",
    use_meta_expression: bool = False,
) -> list[Row]:
    """
    Build source feature rows for DDBJ MSS annotation (wgs_maker).

    Returns a list of 5-column rows.

    When *use_meta_expression* is True the feature is built for inclusion in the
    COMMON block: location is ``1..E``, ``submitter_seqid`` is ``@@[entry]@@``,
    and ``ff_definition`` uses ``@@[...]@@`` meta-notation.  *source_modifier_key*
    (from ``INFRASPECIFIC_NAME_MODIFIER``, e.g. ``"strain"``, ``"cultivar"``) is used
    to select the modifier placeholder in ``ff_definition``; when empty only
    ``@@[organism]@@`` is included.
    """
    if use_meta_expression:
        return _create_source_with_meta(
            _trad_submission_category, source_dict, source_modifier_key
        )

    submitter_seqid = None
    environmental_sample = False
    plasmid = False
    if _trad_submission_category == "GNM":
        mol_type = "genomic DNA"
        if seq_type in ["c", "complete"]:
            ff_def = "@@[organism]@@ @@[strain]@@ DNA, complete genome"
        elif seq_type in ["n", "nearly complete", "nearly-complete"]:
            ff_def = "@@[organism]@@ @@[strain]@@ DNA, nearly complete genome"
        elif seq_type in ["p", "plasmid"]:
            ff_def = "@@[organism]@@ @@[strain]@@ plasmid @@[plasmid]@@ DNA, complete sequence"
            plasmid = True
        else:
            submitter_seqid = "@@[entry]@@"
            ff_def = "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"
    elif _trad_submission_category == "MAG":
        mol_type = "genomic DNA"
        environmental_sample = True
        if seq_type in ["c", "complete"]:
            ff_def = "@@[organism]@@ @@[isolate]@@ DNA, complete genome"
        elif seq_type in ["n", "nearly complete", "nearly-complete"]:
            ff_def = "@@[organism]@@ @@[isolate]@@ DNA, nearly complete genome"
        elif seq_type in ["p", "plasmid"]:
            ff_def = "@@[organism]@@ @@[isolate]@@ plasmid @@[plasmid]@@ DNA, complete sequence"
            plasmid = True
        else:
            submitter_seqid = "@@[entry]@@"
            ff_def = "@@[organism]@@ @@[isolate]@@ DNA, @@[submitter_seqid]@@"
    elif _trad_submission_category == "WGS":
        mol_type = "genomic DNA"
        submitter_seqid = "@@[entry]@@"
        ff_def = "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"
    elif _trad_submission_category == "MAG-WGS":
        mol_type = "genomic DNA"
        environmental_sample = True
        submitter_seqid = "@@[entry]@@"
        ff_def = "@@[organism]@@ @@[isolate]@@ DNA, @@[submitter_seqid]@@"

    ret: list[Row] = []
    ret.append(["", "source", "1..E", "mol_type", mol_type])
    ret.append(["", "", "", "ff_definition", ff_def])
    if submitter_seqid:
        ret.append(["", "", "", "submitter_seqid", submitter_seqid])
    if environmental_sample:
        ret.append(["", "", "", "environmental_sample", ""])
    if plasmid:
        ret.append(["", "", "", "plasmid", seq_name])
    for key, value in source_dict.items():
        ret.append(["", "", "", key, value])
    if _trad_submission_category in ["WGS", "MAG-WGS"]:
        # The source feature will be appended to COMMON. Nothing to do.
        pass
    elif _trad_submission_category in ["GNM", "MAG"]:
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
    environmental_sample = category in ("MAG", "MAG-WGS")
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
    for key, value in source_dict.items():
        if key != "mol_type":
            rows.append(["", "", "", key, str(value)])
    return rows
