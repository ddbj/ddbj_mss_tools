"""
Write DDBJ MSS annotation (.ann) files.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from common.common_builder import create_common
from common.gap_annotator import GapAnnotator, annotate_gaps
from common.fasta import parse_fasta_sequences
from common.source_builder import (
    ChromosomeEntry,
    load_chromosomes,
    source_qualifier,
    ff_definition,
)
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

def _common_placeholder() -> list[Row]:
    """Return placeholder COMMON rows when no common JSON is provided."""
    rows: list[Row] = [
        ["COMMON", "DBLINK",     "", "project",    ""],
        ["",       "",           "", "biosample",   ""],
        ["",       "SUBMITTER",  "", "ab_name",     ""],
        ["",       "",           "", "contact",     ""],
        ["",       "",           "", "email",       ""],
        ["",       "",           "", "institute",   ""],
        ["",       "",           "", "country",     ""],
        ["",       "REFERENCE",  "", "title",       ""],
        ["",       "",           "", "ab_name",     ""],
        ["",       "",           "", "status",      "Unpublished"],
        ["",       "",           "", "year",        ""],
        ["",       "DATE",       "", "hold_date",   ""],
    ]
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
    sequences = parse_fasta_sequences(fsa_path)
    lengths = {seq_id: len(seq) for seq_id, seq in sequences.items()}

    gap_cfg = common.ASSEMBLY_GAP if common is not None else None
    gap_annotators: list[GapAnnotator] = []
    if gap_cfg:
        cfgs = gap_cfg if isinstance(gap_cfg, list) else [gap_cfg]
        gap_annotators = [
            GapAnnotator(
                linkage_evidence=cfg.linkage_evidence,
                min_gap_length=cfg.min_gap_length,
                max_gap_length=cfg.max_gap_length,
                gap_type=cfg.gap_type,
                estimated_length=cfg.estimated_length,
            )
            for cfg in cfgs
            if cfg.enabled
        ]
    if not gap_annotators:
        sequences = {}

    # Use FASTA entry order as the canonical list; fall back to .tbl order for
    # any entries that appear only in the feature table.
    all_ids = list(lengths.keys()) + [e for e in entries if e not in lengths]

    # Base source qualifiers come entirely from common.SOURCE.
    base_source: dict[str, str] = {}
    if common is not None and common.SOURCE:
        base_source.update(common.SOURCE)

    # infraspecific_name_modifier: value of the qualifier named by INFRASPECIFIC_NAME_MODIFIER
    organism = base_source.get("organism", "")
    source_id_key = common.INFRASPECIFIC_NAME_MODIFIER if common is not None else None
    infraspecific_name_modifier = base_source.get(source_id_key, "") if source_id_key else ""

    # WGS mode: all entries are unplaced (not listed in chromosomes, or all type==unplaced)
    def _is_unplaced(eid: str) -> bool:
        if chromosomes is None:
            return True
        e = chromosomes.get(eid)
        return e is None or e.type == "unplaced"

    is_wgs = all(_is_unplaced(eid) for eid in all_ids)

    rows: list[Row] = []
    if common is None:
        rows.extend(_common_placeholder())
    else:
        rows.extend(create_common(common.model_dump(exclude_none=True)))

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
        source_quals.update(source_qualifier(chr_entry, entry_id, is_wgs))
        source_quals["ff_definition"] = ff_definition(
            chr_entry, entry_id, organism, infraspecific_name_modifier, is_wgs
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

        if gap_annotators and entry_id in sequences:
            rows.extend(annotate_gaps(gap_annotators, sequences[entry_id]))

    with open(ann_path, "w") as fout:
        for row in rows:
            fout.write("\t".join(row) + "\n")

    print(f"[convert] → {ann_path}", file=sys.stderr)
