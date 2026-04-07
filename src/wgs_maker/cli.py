"""
wgs_maker — Convert FASTA files to DDBJ MSS submission format.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from .core import create_mss
from .gap_annotator import GapAnnotator
from .schema_util import get_local_schema, load_json_file

LINKAGE_EVIDENCES = [
    "pcr", "paired-ends", "align_genus", "align_xgenus", "align_trnscpt",
    "within_clone", "clone_contig", "map", "strobe", "proximity_ligation",
    "unspecified",
]
GAP_TYPES = [
    "auto", "between_scaffolds", "within_scaffold", "telomere", "centromere",
    "short_arm", "heterochromatin", "repeat_within_scaffold",
    "repeat between_scaffolds", "contamination", "unknown",
]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="wgs_maker",
        description="Convert FASTA file to MSS format for DDBJ submission",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--excel", type=str, help="Excel workbook file")
    group.add_argument("--tsv",   type=str, help="TSV file")

    parser.add_argument(
        "--sheet", type=str, default="Sheet1",
        help="Sheet name in the Excel workbook (default: Sheet1)",
    )
    parser.add_argument(
        "-m", "--metadata_json_file",
        help="Common metadata in JSON format (for submitter and reference)",
    )
    parser.add_argument("-o", "--out_dir", default=".", help="Output directory (default: .)")
    parser.add_argument("-H", "--hold_date", help='Hold date, format "yyyymmdd"')
    parser.add_argument(
        "--linkage_evidence", choices=LINKAGE_EVIDENCES, default="paired-ends",
        help='Linkage evidence for assembly_gap features (default: paired-ends)',
    )
    parser.add_argument(
        "--gap_type", choices=GAP_TYPES, default="auto",
        help="Gap type for assembly_gap features (default: auto)",
    )
    parser.add_argument(
        "--gap_length", choices=["auto", "known", "unknown"], default="auto",
        help="Estimated gap length (default: auto)",
    )
    parser.add_argument(
        "--min_gap_length", type=int, default=10,
        help="Minimum gap length (default: 10)",
    )
    return parser


def main() -> None:
    parser = _build_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    gap_annotator = GapAnnotator.initialize(args)
    base_json_data = load_json_file(args.metadata_json_file)
    base_schema    = get_local_schema()

    if args.excel:
        df = pd.read_excel(args.excel, sheet_name=args.sheet, header=[0, 1], dtype=str)
    else:
        df = pd.read_csv(args.tsv, sep="\t", header=[0, 1], dtype=str)
    df.fillna("", inplace=True)

    for _, row in df.iterrows():
        create_mss(
            row, base_json_data, base_schema, args.out_dir,
            gap_annotator, hold_date=args.hold_date,
        )
