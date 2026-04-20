"""
egapx2mss — Convert NCBI ASN.1 (EGAPx output) to DDBJ MSS submission format.

Workflow
--------
1. Run asn2gb -f t  → NCBI feature table (.tbl)
2. Run asn2fsa      → FASTA sequence file (.fa)
3. Parse both files → DDBJ MSS annotation (.ann) + cleaned FASTA (.fa)

The tools asn2gb and asn2fsa are downloaded automatically from
https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/cmdline/ if they are not found
in --bin-dir (default: ~/.local/share/ddbj_mss_tools/bin).
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Optional

from .ann_writer import write_ddbj_ann
from common.source_builder import load_chromosomes
from .asn_tools import DEFAULT_BIN_DIR, ensure_tools, run_asn2fsa, run_asn2gb_tbl
from common.fasta import write_clean_fasta
from .models import CommonModel, load_common_json


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="egapx2mss",
        description="Convert NCBI ASN.1 (EGAPx) to DDBJ MSS format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input", help="Input ASN.1 file (.asn)")
    parser.add_argument(
        "-o", "--output",
        help="Output file prefix (default: input basename without extension)",
    )
    parser.add_argument(
        "--bin-dir", default=None,
        help=(
            "Directory that contains (or should receive) asn2gb and asn2fsa binaries. "
            f"Default: {DEFAULT_BIN_DIR}"
        ),
    )
    parser.add_argument(
        "--common",
        help="JSON file with common submission metadata (DBLINK, SUBMITTER, REFERENCE, DATE)",
    )
    parser.add_argument(
        "--chromosomes",
        help=(
            "3-column TSV file mapping entry names to source qualifiers "
            "(entry_name <TAB> qualifier_key <TAB> qualifier_value). "
            "Entries absent from this file get submitter_seqid set to the entry name."
        ),
    )
    parser.add_argument(
        "--keep-tmp", action="store_true",
        help="Keep intermediate .tbl and raw FASTA files",
    )
    args = parser.parse_args()

    # ── Resolve paths ──────────────────────────────────────────────────────
    asn_path = os.path.abspath(args.input)
    bin_dir  = Path(args.bin_dir).resolve() if args.bin_dir else DEFAULT_BIN_DIR

    prefix      = args.output if args.output else os.path.splitext(asn_path)[0]
    out_tbl     = prefix + ".tbl"
    out_fsa_raw = prefix + "_raw.fa"
    out_fsa     = prefix + ".fa"
    out_ann     = prefix + ".ann"

    # ── Ensure tools ──────────────────────────────────────────────────────
    asn2gb, asn2fsa = ensure_tools(bin_dir)

    # ── Step 1: feature table ─────────────────────────────────────────────
    print("[step 1/3] Generating feature table (asn2gb) ...", file=sys.stderr)
    run_asn2gb_tbl(asn2gb, asn_path, out_tbl)

    # ── Step 2: FASTA ────────────────────────────────────────────────────
    print("[step 2/3] Generating FASTA (asn2fsa) ...", file=sys.stderr)
    tmpdir = tempfile.mkdtemp(prefix="egapx2mss_")
    try:
        run_asn2fsa(asn2fsa, asn_path, out_fsa_raw, tmpdir)
    finally:
        if not args.keep_tmp:
            shutil.rmtree(tmpdir, ignore_errors=True)

    write_clean_fasta(out_fsa_raw, out_fsa)
    if not args.keep_tmp:
        os.unlink(out_fsa_raw)

    # ── Step 3: DDBJ MSS annotation ──────────────────────────────────────
    print("[step 3/3] Converting to DDBJ MSS annotation ...", file=sys.stderr)
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

    write_ddbj_ann(
        out_tbl, out_fsa, out_ann,
        common=common,
        chromosomes=chromosomes,
    )

    if not args.keep_tmp:
        os.unlink(out_tbl)

    print("\nDone.", file=sys.stderr)
    print(f"  Annotation : {out_ann}", file=sys.stderr)
    print(f"  FASTA      : {out_fsa}", file=sys.stderr)
