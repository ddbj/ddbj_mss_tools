"""
Write DDBJ MSS annotation (.ann) files from a FASTA file only (no feature table).

Used by mss_builder. For assemblies that only need source + assembly_gap features.
"""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING, Optional

from common.common_builder import create_common
from common.fasta import parse_fasta_sequences
from common.gap_annotator import GapAnnotator
from common.source_builder import (
    ChromosomeEntry,
    ff_definition,
    source_qualifier,
)

if TYPE_CHECKING:
    from common.models import CommonModel

Row = list[str]


# ── Placeholder COMMON sections ───────────────────────────────────────────────

def _common_placeholder_wgs() -> list[Row]:
    """Placeholder COMMON rows for WGS mode (source included in COMMON)."""
    return [
        ["COMMON", "DBLINK",    "", "project",           ""],
        ["",       "",          "", "biosample",          ""],
        ["",       "SUBMITTER", "", "ab_name",            ""],
        ["",       "",          "", "contact",            ""],
        ["",       "",          "", "email",              ""],
        ["",       "",          "", "institute",          ""],
        ["",       "",          "", "country",            ""],
        ["",       "REFERENCE", "", "title",              ""],
        ["",       "",          "", "ab_name",            ""],
        ["",       "",          "", "status",             "Unpublished"],
        ["",       "",          "", "year",               ""],
        ["",       "DATE",      "", "hold_date",          ""],
        ["",       "source",    "1..E", "mol_type",       "genomic DNA"],
        ["",       "",          "", "ff_definition",      "@@[organism]@@ DNA, @@[submitter_seqid]@@"],
        ["",       "",          "", "submitter_seqid",    "@@[entry]@@"],
        ["",       "",          "", "organism",           ""],
    ]


def _common_placeholder_nonwgs() -> list[Row]:
    """Placeholder COMMON rows for non-WGS mode (source written per entry)."""
    return [
        ["COMMON", "DBLINK",    "", "project",     ""],
        ["",       "",          "", "biosample",    ""],
        ["",       "SUBMITTER", "", "ab_name",      ""],
        ["",       "",          "", "contact",      ""],
        ["",       "",          "", "email",        ""],
        ["",       "",          "", "institute",    ""],
        ["",       "",          "", "country",      ""],
        ["",       "REFERENCE", "", "title",        ""],
        ["",       "",          "", "ab_name",      ""],
        ["",       "",          "", "status",       "Unpublished"],
        ["",       "",          "", "year",         ""],
        ["",       "DATE",      "", "hold_date",    ""],
    ]


# ── Main writer ───────────────────────────────────────────────────────────────

def write_mss_ann(
    fsa_path: str,
    ann_path: str,
    common: Optional["CommonModel"] = None,
    chromosomes: Optional[dict[str, ChromosomeEntry]] = None,
) -> None:
    """
    Parse *fsa_path* (FASTA) and write a DDBJ MSS annotation file to *ann_path*.

    If *common* is provided (a validated CommonModel), its values are written
    into the COMMON section; otherwise placeholder lines are written.

    For WGS submissions (no *chromosomes* provided, or all entries unplaced):
    - The source feature is placed in the COMMON block using ``@@[entry]@@`` and
      ``@@[organism]@@ DNA, @@[submitter_seqid]@@`` meta-notation.
    - Per-entry body contains only assembly_gap features (if any).

    For non-WGS submissions (*chromosomes* provided with placed sequences):
    - Source features are written per entry with chromosome/organelle qualifiers.
    - Assembly_gap features follow the source for each entry.
    """
    sequences = parse_fasta_sequences(fsa_path)
    lengths = {seq_id: len(seq) for seq_id, seq in sequences.items()}
    all_ids = list(lengths.keys())

    gap_cfg = common.ASSEMBLY_GAP if common is not None else None
    if not gap_cfg:
        sequences = {}
    gap_annotator = (
        GapAnnotator(
            linkage_evidence=gap_cfg.linkage_evidence,
            min_gap_length=gap_cfg.min_gap_length,
        )
        if gap_cfg else None
    )

    # Determine WGS mode: no chromosomes file, or every entry is unplaced
    def _is_unplaced(eid: str) -> bool:
        if chromosomes is None:
            return True
        e = chromosomes.get(eid)
        return e is None or e.type == "unplaced"

    is_wgs = all(_is_unplaced(eid) for eid in all_ids)

    # Base source qualifiers from common.SOURCE
    base_source: dict[str, str] = {}
    if common is not None and common.SOURCE:
        base_source.update(common.SOURCE)

    organism = base_source.get("organism", "")
    source_id_key = common.SOURCE_MODIFIER if common is not None else None
    source_modifier = base_source.get(source_id_key, "") if source_id_key else ""

    rows: list[Row] = []

    # ── COMMON section ────────────────────────────────────────────────────────
    if common is None:
        if is_wgs:
            rows.extend(_common_placeholder_wgs())
        else:
            rows.extend(_common_placeholder_nonwgs())
    else:
        common_dict = common.model_dump(exclude_none=True)
        if is_wgs:
            # Inject category so _build_common_source picks the right template
            common_dict["_trad_submission_category"] = "WGS"
            rows.extend(create_common(common_dict, include_source=True))
        else:
            rows.extend(create_common(common_dict))

    # ── Per-entry body ────────────────────────────────────────────────────────
    for entry_id in all_ids:
        length = lengths[entry_id]
        location = f"1..{length}"

        chr_entry: Optional[ChromosomeEntry] = (
            chromosomes.get(entry_id) if chromosomes else None
        )
        is_circular = chr_entry.is_circular if chr_entry is not None else False

        if is_circular:
            rows.append([entry_id, "TOPOLOGY", "", "circular", ""])

        if not is_wgs:
            # Source feature per entry
            source_quals: dict[str, str] = dict(base_source)
            source_quals.update(source_qualifier(chr_entry, entry_id, is_wgs=False))
            source_quals["ff_definition"] = ff_definition(
                chr_entry, entry_id, organism, source_modifier, is_wgs=False
            )

            source_entry_col = "" if is_circular else entry_id
            qual_items = list(source_quals.items())
            first_key, first_val = qual_items[0]
            rows.append([source_entry_col, "source", location, first_key, first_val])
            for q_key, q_val in qual_items[1:]:
                rows.append(["", "", "", q_key, q_val])

        # Assembly gap features
        if gap_annotator and entry_id in sequences:
            # In WGS mode the first gap row carries the entry name (since source
            # is in COMMON). In non-WGS mode source already opened the entry, so
            # gap rows use an empty first column.
            seq_name_for_gap = entry_id if is_wgs and not is_circular else None
            gap_rows = gap_annotator.annotate(sequences[entry_id], seq_name=seq_name_for_gap)
            rows.extend(gap_rows)

    with open(ann_path, "w") as fout:
        for row in rows:
            fout.write("\t".join(row) + "\n")

    print(f"[mss_builder] → {ann_path}", file=sys.stderr)
