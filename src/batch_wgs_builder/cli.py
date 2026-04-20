"""
batch_wgs_builder — Batch-create DDBJ MSS WGS/MAG-WGS submission files.
"""

from __future__ import annotations

import argparse
import sys

from .core import run


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="batch_wgs_builder",
        description="Batch-create DDBJ MSS WGS/MAG-WGS submission files from a sample list TSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "tsv",
        help="Sample list TSV (2-row header: row1=feature names, row2=qualifier names)",
    )
    parser.add_argument(
        "--common", "-m",
        help="Common metadata JSON (SUBMITTER, REFERENCE, ASSEMBLY_GAP, SOURCE_MODIFIER, …)",
    )
    parser.add_argument(
        "--out-dir", "-o", default=".",
        help="Output directory (default: current directory)",
    )
    parser.add_argument(
        "--hold-date", "-H",
        help='Hold date in YYYYMMDD format (written as DATE.hold_date in annotation)',
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    try:
        run(
            tsv_path=args.tsv,
            common_path=args.common,
            out_dir=args.out_dir,
            hold_date=args.hold_date,
        )
    except Exception as exc:
        sys.exit(f"Error: {exc}")

    print("\nDone.", file=sys.stderr)
