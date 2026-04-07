"""
Write DDBJ MSS annotation (.ann) files.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING, Optional

from .fasta import parse_fasta_lengths, parse_fasta_sequences
from .tbl_parser import collect_qualifiers, format_location, parse_tbl

if TYPE_CHECKING:
    from .models import CommonModel

Row = list[str]


# ── Internal helpers ──────────────────────────────────────────────────────────

def _lookup_length(entry_id: str, lengths: dict[str, int]) -> int:
    """Return sequence length for *entry_id*, with fallback to partial-match."""
    if entry_id in lengths:
        return lengths[entry_id]
    for sid, slen in lengths.items():
        if entry_id in sid or sid in entry_id:
            return slen
    return 0


def _feature_rows(first_col: str, feat: dict) -> list[Row]:
    """Return rows for one feature block."""
    quals = collect_qualifiers(feat)
    if not quals:
        return []

    ftype    = feat["type"]
    location = format_location(feat["intervals"])

    rows: list[Row] = []
    q0_name, q0_val = quals[0]
    rows.append([first_col, ftype, location, q0_name, q0_val if q0_val is not None else ""])

    for q_name, q_val in quals[1:]:
        rows.append(["", "", "", q_name, q_val if q_val is not None else ""])

    return rows


# ── COMMON section builder ────────────────────────────────────────────────────

def build_common(common: Optional["CommonModel"]) -> list[Row]:
    """Return rows for the COMMON block."""
    rows: list[Row] = []

    if common is not None and common.DATATYPE:
        items = list(common.DATATYPE.items())
        first_key, first_val = items[0]
        rows.append(["", "DATATYPE", "", first_key, first_val])
        for key, val in items[1:]:
            rows.append(["", "", "", key, val])

    if common is not None and common.KEYWORD:
        first_row = True
        for key, val in common.KEYWORD.items():
            values = [val] if isinstance(val, str) else val
            for v in values:
                if first_row:
                    rows.append(["", "KEYWORD", "", key, v])
                    first_row = False
                else:
                    rows.append(["", "", "", key, v])

    if common is None:
        rows.append(["", "DBLINK", "", "project", ""])
        rows.append(["", "", "", "biosample", ""])
        rows.append(["", "SUBMITTER", "", "ab_name", ""])
        rows.append(["", "", "", "contact", ""])
        rows.append(["", "", "", "email", ""])
        rows.append(["", "", "", "institute", ""])
        rows.append(["", "", "", "country", ""])
        rows.append(["", "REFERENCE", "", "title", ""])
        rows.append(["", "", "", "ab_name", ""])
        rows.append(["", "", "", "status", "Unpublished"])
        rows.append(["", "", "", "year", ""])
        rows.append(["", "DATE", "", "hold_date", ""])
        rows[0][0] = "COMMON"
        return rows

    # DBLINK
    db = common.DBLINK
    rows.append(["", "DBLINK", "", "project", db.project])
    rows.append(["", "", "", "biosample", db.sample])
    if db.DRA:
        for acc in db.DRA:
            rows.append(["", "", "", "DRA", acc])

    # SUBMITTER
    sub = common.SUBMITTER
    if sub:
        first_name, *rest_names = sub.ab_name
        rows.append(["", "SUBMITTER", "", "ab_name", first_name])
        for name in rest_names:
            rows.append(["", "", "", "ab_name", name])
        for field in ("contact", "email", "phone", "fax", "institute",
                      "department", "country", "state", "city", "street", "zip"):
            val = getattr(sub, field)
            if val is not None:
                rows.append(["", "", "", field, val])
    else:
        rows.append(["", "SUBMITTER", "", "ab_name", ""])

    # REFERENCE
    refs = common.REFERENCE
    if refs:
        for ref in refs:
            first_name, *rest_names = ref.ab_name
            rows.append(["", "REFERENCE", "", "title", ref.title])
            rows.append(["", "", "", "ab_name", first_name])
            for name in rest_names:
                rows.append(["", "", "", "ab_name", name])
            rows.append(["", "", "", "status", ref.status])
            rows.append(["", "", "", "year", str(ref.year)])
            for field in ("journal", "volume", "start_page", "end_page"):
                val = getattr(ref, field)
                if val is not None:
                    rows.append(["", "", "", field, val])
    else:
        rows.append(["", "REFERENCE", "", "title", ""])
        rows.append(["", "", "", "ab_name", ""])
        rows.append(["", "", "", "status", "Unpublished"])
        rows.append(["", "", "", "year", ""])

    # DATE
    if common.DATE:
        hold = common.DATE.get("hold_date", "")
        rows.append(["", "DATE", "", "hold_date", hold])
    else:
        rows.append(["", "DATE", "", "hold_date", ""])

    if rows:
        rows[0][0] = "COMMON"

    return rows


# ── Chromosome table ─────────────────────────────────────────────────────────

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


def _source_qualifier(entry: Optional[ChromosomeEntry], seq_id: str,
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


def _ff_definition(entry: Optional[ChromosomeEntry], seq_id: str, organism: str,
                   source_identifier: str, is_wgs: bool = False) -> str:
    """
    Build the ff_definition qualifier value following mss_format.md.

    *source_identifier* is the value of the qualifier named by SOURCE_IDENTIFIER
    (e.g. the value of 'strain' or 'isolate') from common.SOURCE.

    *is_wgs* is True when all entries in the submission are unplaced (WGS mode).
    """
    prefix = f"{organism} {source_identifier}".strip() if source_identifier else organism

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


# ── Assembly gap detection ────────────────────────────────────────────────────

# gap_type and gap_length are determined by linkage_evidence
_GAP_ATTRS: dict[str, tuple[str, str]] = {
    "paired-ends":        ("within scaffolds", "known"),
    "proximity ligation": ("within scaffolds", "unknown"),
    "align genus":        ("within scaffolds", "unknown"),
}


def _gap_rows(seq: str, linkage_evidence: str, min_gap_length: int) -> list[Row]:
    """Return assembly_gap feature rows for all N-runs in *seq*."""
    gap_type, estimated_length = _GAP_ATTRS[linkage_evidence]
    pattern = re.compile(f"[Nn]{{{min_gap_length},}}")
    rows: list[Row] = []
    for match in pattern.finditer(seq):
        start = match.start() + 1  # 1-based
        end = match.end()
        rows.append(["", "assembly_gap", f"{start}..{end}", "estimated_length", estimated_length])
        rows.append(["", "", "", "gap_type", gap_type])
        rows.append(["", "", "", "linkage_evidence", linkage_evidence])
    return rows


# ── Main writer ───────────────────────────────────────────────────────────────

def write_ddbj_ann(
    tbl_path: str,
    fsa_path: str,
    ann_path: str,
    common: Optional["CommonModel"] = None,
    chromosomes: Optional[dict[str, ChromosomeEntry]] = None,
) -> None:
    """
    Parse *tbl_path* (NCBI feature table) and *fsa_path* (FASTA),
    then write a DDBJ MSS annotation file to *ann_path*.

    If *common* is provided (a validated CommonModel), its values are written
    into the COMMON section; otherwise placeholder lines are written.
    Source feature qualifiers are taken from ``common.SOURCE`` if present.
    Assembly gap detection is driven by ``common.ASSEMBLY_GAP`` if present.

    If *chromosomes* is provided (entry_name → (qualifier_key, qualifier_value, is_circular)),
    the matching qualifier overrides the corresponding key in the source feature.
    Entries absent from the table get ``submitter_seqid`` set to the entry name.
    """
    import sys

    entries = parse_tbl(tbl_path)
    lengths = parse_fasta_lengths(fsa_path)

    gap_cfg = common.ASSEMBLY_GAP if common is not None else None
    sequences = parse_fasta_sequences(fsa_path) if gap_cfg else {}

    # Use FASTA entry order as the canonical list; fall back to .tbl order for
    # any entries that appear only in the feature table.
    all_ids = list(lengths.keys()) + [e for e in entries if e not in lengths]

    # Base source qualifiers come entirely from common.SOURCE.
    base_source: dict[str, str] = {}
    if common is not None and common.SOURCE:
        base_source.update(common.SOURCE)

    # source_identifier: value of the qualifier named by SOURCE_IDENTIFIER
    organism = base_source.get("organism", "")
    source_id_key = common.SOURCE_IDENTIFIER if common is not None else None
    source_identifier = base_source.get(source_id_key, "") if source_id_key else ""

    # WGS mode: all entries are unplaced (not listed in chromosomes, or all type==unplaced)
    def _is_unplaced(eid: str) -> bool:
        if chromosomes is None:
            return True
        e = chromosomes.get(eid)
        return e is None or e.type == "unplaced"

    is_wgs = all(_is_unplaced(eid) for eid in all_ids)

    rows: list[Row] = []
    rows.extend(build_common(common))

    for entry_id in all_ids:
        length = _lookup_length(entry_id, lengths)
        location = f"1..{length}" if length else "1.."

        # Resolve ChromosomeEntry (None = unplaced/default)
        chr_entry: Optional[ChromosomeEntry] = chromosomes.get(entry_id) if chromosomes else None

        is_circular = chr_entry.is_circular if chr_entry is not None else False

        if is_circular:
            rows.append([entry_id, "TOPOLOGY", "", "circular", ""])

        # Build per-entry source qualifiers: base + entry-specific + ff_definition
        source_quals: dict[str, str] = dict(base_source)
        source_quals.update(_source_qualifier(chr_entry, entry_id, is_wgs))
        source_quals["ff_definition"] = _ff_definition(
            chr_entry, entry_id, organism, source_identifier, is_wgs
        )

        # source feature: entry_id on the TOPOLOGY row if circular, else on source row
        source_entry_col = "" if is_circular else entry_id
        qual_items = list(source_quals.items())
        first_key, first_val = qual_items[0]
        rows.append([source_entry_col, "source", location, first_key, first_val])
        for q_key, q_val in qual_items[1:]:
            rows.append(["", "", "", q_key, q_val])

        for feat in entries.get(entry_id, []):
            rows.extend(_feature_rows("", feat))

        if gap_cfg and entry_id in sequences:
            rows.extend(_gap_rows(sequences[entry_id], gap_cfg.linkage_evidence, gap_cfg.min_gap_length))

    with open(ann_path, "w") as fout:
        for row in rows:
            fout.write("\t".join(row) + "\n")

    print(f"[convert] → {ann_path}", file=sys.stderr)
