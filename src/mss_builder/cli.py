"""
mss_builder — Create DDBJ MSS submission files from a FASTA file.

Workflow
--------
1. Read input FASTA; write clean copy with '//' separators → .fa
2. Build DDBJ MSS annotation → .ann
   - COMMON block from --common JSON (or placeholder if omitted)
   - WGS mode (default, no --chromosomes): source in COMMON with @@[entry]@@ / @@[organism]@@ meta-notation
   - Chromosome mode (--chromosomes provided): per-entry source with chromosome/organelle qualifiers
   - assembly_gap features detected from N-runs (requires ASSEMBLY_GAP in --common JSON)
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Optional

from common.fasta import write_clean_fasta
from common.models import CommonModel, load_common_json
from common.source_builder import load_chromosomes

from .ann_writer import write_mss_ann


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="mss_builder",
        description="Create DDBJ MSS submission files from a FASTA file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input", help="Input FASTA file (.fa / .fasta)")
    parser.add_argument(
        "-o", "--output",
        help="Output file prefix (default: input basename without extension)",
    )
    parser.add_argument(
        "--common",
        help=(
            "JSON file with common submission metadata "
            "(DBLINK, SUBMITTER, REFERENCE, DATE, SOURCE, ASSEMBLY_GAP, …). "
            "If omitted, placeholder lines are written."
        ),
    )
    parser.add_argument(
        "--chromosomes",
        help=(
            "5-column TSV: seq_id <TAB> type <TAB> seq_name <TAB> status <TAB> topology. "
            "type is one of: chromosome, organelle, unplaced. "
            "If omitted, all sequences are treated as WGS (unplaced) contigs."
        ),
    )
    args = parser.parse_args()

    # ── Resolve paths ─────────────────────────────────────────────────────────
    fasta_path = os.path.abspath(args.input)
    prefix = args.output if args.output else os.path.splitext(fasta_path)[0]
    out_fsa = os.path.abspath(prefix + ".fa")
    out_ann = os.path.abspath(prefix + ".ann")

    # ── Load optional metadata ────────────────────────────────────────────────
    common: Optional[CommonModel] = None
    if args.common:
        try:
            common = load_common_json(args.common)
        except Exception as exc:
            sys.exit(f"Error: failed to load --common file '{args.common}':\n{exc}")

    chromosomes = None
    if args.chromosomes:
        try:
            chromosomes = load_chromosomes(args.chromosomes)
        except Exception as exc:
            sys.exit(f"Error: failed to load --chromosomes file '{args.chromosomes}':\n{exc}")

    # ── Step 1: clean FASTA ───────────────────────────────────────────────────
    print("[step 1/2] Writing clean FASTA with '//' separators ...", file=sys.stderr)
    if os.path.abspath(fasta_path) == out_fsa:
        # Input and output are the same path — use a temp file to avoid clobbering
        tmp_fd, tmp_path = tempfile.mkstemp(suffix=".fa", prefix="mss_builder_")
        os.close(tmp_fd)
        try:
            write_clean_fasta(fasta_path, tmp_path)
            shutil.move(tmp_path, out_fsa)
        except Exception:
            os.unlink(tmp_path)
            raise
    else:
        write_clean_fasta(fasta_path, out_fsa)

    # ── Step 2: DDBJ MSS annotation ───────────────────────────────────────────
    print("[step 2/2] Building DDBJ MSS annotation ...", file=sys.stderr)
    write_mss_ann(out_fsa, out_ann, common=common, chromosomes=chromosomes)

    print("\nDone.", file=sys.stderr)
    print(f"  Annotation : {out_ann}", file=sys.stderr)
    print(f"  FASTA      : {out_fsa}", file=sys.stderr)
