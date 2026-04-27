"""
batch_wgs_builder — Batch-create DDBJ MSS WGS/MAG-WGS submission files.

Workflow
--------
1. Read common JSON (SUBMITTER, REFERENCE, ASSEMBLY_GAP, INFRASPECIFIC_NAME_MODIFIER, …)
2. Read sample_list TSV (2-row header: row0=feature names, row1=qualifier names)
3. For each sample row:
   - Merge DBLINK, ST_COMMENT, source, COMMENT from TSV with common JSON
   - Inject DATATYPE/KEYWORD/DIVISION based on submission category
   - Write COMMON block (source feature in COMMON with @@[entry]@@ meta-notation)
   - Add per-sequence assembly_gap features
   - Output {prefix}.ann and {prefix}.fa
"""

from __future__ import annotations

import copy
import os
import re
import sys

import pandas as pd

from common.common_builder import create_common
from common.fasta import read_fasta
from common.gap_annotator import GapAnnotator, annotate_gaps
from common.models import AssemblyGapModel

Row = list[str]

_KEYWORDS_WGS = ["WGS", "STANDARD_DRAFT"]
_KEYWORDS_MAG_WGS = ["ENV", "WGS", "STANDARD_DRAFT", "Metagenome Assembled Genome", "MAG"]


def parse_tsv(tsv_path: str) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t", header=[0, 1], dtype=str)
    df.fillna("", inplace=True)
    return df


def _split_values(value: str) -> list[str]:
    """Split semicolon- or comma-separated string into a list of stripped values."""
    return [v.strip() for v in re.split(r"[;,]", value) if v.strip()]


def build_sample_json(row: pd.Series, common_base: dict) -> tuple[str, str, dict]:
    """
    Build (file_path, category, common_json) from one TSV row merged with common_base.

    common_base comes from the common JSON and may contain SUBMITTER, REFERENCE,
    SOURCE (as defaults), INFRASPECIFIC_NAME_MODIFIER, ASSEMBLY_GAP, etc.
    TSV row data overrides common_base values where they overlap.
    """
    file_path = ""
    category = "WGS"
    dblink: dict = copy.deepcopy(common_base.get("DBLINK", {}))
    st_comment: dict = copy.deepcopy(common_base.get("ST_COMMENT", {}))
    source_dict: dict = copy.deepcopy(common_base.get("SOURCE", {}))
    comment_lines: list[str] = []

    for (feature, qualifier) in row.index:
        value = str(row[(feature, qualifier)]).strip()
        if not value:
            continue

        if feature == "_":
            if qualifier == "_file_path":
                file_path = value
            elif qualifier == "_trad_submission_category":
                category = value
        elif feature == "_sequence":
            pass  # seq_names/types/topologies not used for WGS/MAG-WGS
        elif feature == "source":
            source_dict[qualifier] = value
        elif feature == "COMMENT":
            lines = [v.strip() for v in value.split(";") if v.strip()]
            comment_lines.extend(lines if lines else [value])
        elif feature == "DBLINK":
            if qualifier in ("biosample", "sequence read archive"):
                dblink[qualifier] = _split_values(value)
            else:
                dblink[qualifier] = value
        elif feature == "ST_COMMENT":
            st_comment[qualifier] = value

    # Auto-inject tagset_id for Genome-Assembly-Data if ST_COMMENT has data from TSV
    if st_comment and "tagset_id" not in st_comment:
        st_comment["tagset_id"] = "Genome-Assembly-Data"

    # Build ordered common_json dict (order determines ann file output order)
    common_json: dict = {}

    common_json["DATATYPE"] = {"type": "WGS"}

    if category == "MAG-WGS":
        common_json["DIVISION"] = {"division": "ENV"}
        common_json["KEYWORD"] = {"keyword": _KEYWORDS_MAG_WGS}
    else:
        common_json["KEYWORD"] = {"keyword": _KEYWORDS_WGS}

    if dblink:
        common_json["DBLINK"] = dblink

    if "SUBMITTER" in common_base:
        common_json["SUBMITTER"] = common_base["SUBMITTER"]
    if "REFERENCE" in common_base:
        common_json["REFERENCE"] = common_base["REFERENCE"]
    if "DATE" in common_base:
        common_json["DATE"] = common_base["DATE"]

    if st_comment:
        common_json["ST_COMMENT"] = st_comment

    if comment_lines:
        common_json["COMMENT"] = [{"line": comment_lines}]

    # Tool-specific config keys (used by create_common / source builder)
    common_json["SOURCE"] = source_dict
    common_json["INFRASPECIFIC_NAME_MODIFIER"] = common_base.get("INFRASPECIFIC_NAME_MODIFIER", "")
    common_json["_trad_submission_category"] = category

    return file_path, category, common_json


def _output_prefix(out_dir: str, common_json: dict) -> str:
    """Derive output file prefix from biosample + strain/isolate."""
    dblink = common_json.get("DBLINK", {})
    biosample = dblink.get("biosample", [])
    if isinstance(biosample, list):
        biosample = biosample[0] if biosample else "NO_BIOSAMPLE"

    source = common_json.get("SOURCE", {})
    identifier = source.get("strain") or source.get("isolate") or "NO_IDENTIFIER"
    prefix = f"{biosample}_{identifier}".replace(" ", "_")
    return os.path.join(out_dir, prefix)


def process_sample(
    file_path: str,
    common_json: dict,
    out_dir: str,
    gap_cfg: AssemblyGapModel | list[AssemblyGapModel] | None,
    hold_date: str | None,
) -> None:
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

    ann_rows: list[Row] = create_common(common_json, include_source=True)

    if hold_date:
        ann_rows.append(["", "DATE", "", "hold_date", hold_date])

    seq_records = read_fasta(file_path)
    for seq_record in seq_records:
        entry_id = seq_record.id
        if gap_annotators:
            ann_rows.extend(annotate_gaps(gap_annotators, str(seq_record.seq), seq_name=entry_id))
        seq_record.name = ""
        seq_record.description = ""

    prefix = _output_prefix(out_dir, common_json)
    os.makedirs(out_dir, exist_ok=True)
    out_ann = prefix + ".ann"
    out_fa = prefix + ".fa"

    with open(out_ann, "w") as f:
        for row in ann_rows:
            f.write("\t".join(row) + "\n")

    with open(out_fa, "w") as f:
        for seq_record in seq_records:
            seq_record.seq = seq_record.seq.lower()
            f.write(seq_record.format("fasta"))
            f.write("//\n")

    print(f"  Ann : {out_ann}", file=sys.stderr)
    print(f"  FA  : {out_fa}", file=sys.stderr)


def _load_common_base(path: str) -> tuple[dict, AssemblyGapModel | list[AssemblyGapModel] | None]:
    """Load common JSON without requiring DBLINK (which is per-sample in the TSV)."""
    import json
    from pathlib import Path
    from common.models import _strip_trailing_commas, AssemblyGapModel

    raw = Path(path).read_text(encoding="utf-8")
    data: dict = json.loads(_strip_trailing_commas(raw))

    gap_cfg: AssemblyGapModel | list[AssemblyGapModel] | None = None
    if "ASSEMBLY_GAP" in data:
        raw_gap = data["ASSEMBLY_GAP"]
        if isinstance(raw_gap, list):
            gap_cfg = [AssemblyGapModel.model_validate(item) for item in raw_gap]
        else:
            gap_cfg = AssemblyGapModel.model_validate(raw_gap)

    return data, gap_cfg


def run(
    tsv_path: str,
    common_path: str | None,
    out_dir: str,
    hold_date: str | None = None,
) -> None:
    common_base: dict = {}
    gap_cfg: AssemblyGapModel | None = None

    if common_path:
        common_base, gap_cfg = _load_common_base(common_path)

    df = parse_tsv(tsv_path)

    for i, (_, row) in enumerate(df.iterrows(), 1):
        file_path, category, common_json = build_sample_json(row, common_base)
        if not file_path:
            print(f"[sample {i}] skipping: no _file_path", file=sys.stderr)
            continue
        print(f"[sample {i}] {category}: {file_path}", file=sys.stderr)
        try:
            process_sample(file_path, common_json, out_dir, gap_cfg, hold_date)
        except Exception as exc:
            print(f"[sample {i}] ERROR: {exc}", file=sys.stderr)
            raise
