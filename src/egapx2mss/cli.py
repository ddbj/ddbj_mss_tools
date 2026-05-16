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

When --tbl and --fsa are both provided, steps 1 and 2 are skipped
and no ASN.1 input file is needed.
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
from common.cli_args import validate_prefix, resolve_output
from common.source_builder import load_sequence_roles
from .asn_tools import DEFAULT_BIN_DIR, ensure_tools, run_asn2fsa, run_asn2gb_tbl
from common.fasta import write_clean_fasta
from common.models import CommonModel, load_common_json


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="egapx2mss",
        description="Convert NCBI ASN.1 (EGAPx) to DDBJ MSS format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "input", nargs="?",
        help="Input ASN.1 file (.asn). Omit when --tbl and --fsa are both provided.",
    )
    parser.add_argument(
        "-o", "--outdir", default=None,
        help=(
            "Output directory (created if absent). "
            "Default: same directory as input .asn or --tbl file."
        ),
    )
    parser.add_argument(
        "-p", "--prefix", default=None,
        help=(
            "Output filename prefix (basename only, no directory separators). "
            "Default: basename of input .asn or --tbl file without extension."
        ),
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
        "--submission_category",
        metavar="CATEGORY",
        help=(
            "Submission category (e.g. WGS, GNM, MAG-WGS). "
            "Overrides _submission_category in --common JSON."
        ),
    )
    parser.add_argument(
        "--sequence_roles", "--chromosomes",
        dest="sequence_roles",
        metavar="TSV",
        help=(
            "Sequence role file (5-column TSV): "
            "seq_id <TAB> type <TAB> seq_name <TAB> status <TAB> topology. "
            "type is one of: chromosome, organelle, unplaced. "
            "Entries absent from this file get submitter_seqid set to the entry name. "
            "(--chromosomes is accepted as a legacy alias.)"
        ),
    )
    parser.add_argument(
        "--tbl", default=None,
        help="Pre-existing NCBI feature table (.tbl); skips step 1/3.",
    )
    parser.add_argument(
        "--fsa", default=None,
        help="Pre-existing FASTA file (.fa/.fsa); skips step 2/3.",
    )
    parser.add_argument(
        "--keep-tmp", action="store_true",
        help="Keep intermediate .tbl and raw FASTA files",
    )
    parser.add_argument(
        "--preconvert-only", action="store_true",
        help="Run step 1/3 and step 2/3 only (generate .tbl and .fa); skip MSS annotation conversion",
    )
    args = parser.parse_args()

    # ── Validate argument combinations ────────────────────────────────────
    if bool(args.tbl) != bool(args.fsa):
        parser.error("--tbl and --fsa must be specified together")
    direct_mode = bool(args.tbl and args.fsa)
    if not args.input and not direct_mode:
        parser.error("input .asn file is required unless both --tbl and --fsa are provided")
    if args.input and (args.tbl or args.fsa):
        parser.error("input .asn cannot be combined with --tbl / --fsa")
    validate_prefix(args, parser)
    if args.preconvert_only and direct_mode:
        print("Warning: --preconvert-only has no effect when --tbl and --fsa are both provided.", file=sys.stderr)

    # ── Resolve paths ──────────────────────────────────────────────────────
    bin_dir  = Path(args.bin_dir).resolve() if args.bin_dir else DEFAULT_BIN_DIR
    ref_path = Path(os.path.abspath(args.tbl if direct_mode else args.input))

    out_prefix, out_dir = resolve_output(args, ref_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    base        = out_dir / out_prefix
    out_tbl     = str(base) + ".tbl"
    out_fsa_raw = str(base) + "_raw.fa"
    out_fsa     = str(base) + ".fa"
    out_ann     = str(base) + ".ann"

    # ── Ensure tools (skip when both files are pre-supplied) ──────────────
    if not direct_mode:
        asn2gb, asn2fsa = ensure_tools(bin_dir)
        asn_path = str(ref_path)

    # ── Step 1: feature table ─────────────────────────────────────────────
    if args.tbl:
        out_tbl = os.path.abspath(args.tbl)
        tbl_preexisted = True
        print(f"[step 1/3] Using provided .tbl: {out_tbl}", file=sys.stderr)
    else:
        tbl_preexisted = os.path.exists(out_tbl)
        if tbl_preexisted:
            print(f"[step 1/3] Skipping — {out_tbl} already exists.", file=sys.stderr)
        else:
            print("[step 1/3] Generating feature table (asn2gb) ...", file=sys.stderr)
            run_asn2gb_tbl(asn2gb, asn_path, out_tbl)

    # ── Step 2: FASTA ────────────────────────────────────────────────────
    if args.fsa:
        fsa_source = os.path.abspath(args.fsa)
        fsa_preexisted = True
        if fsa_source == out_fsa:
            print(f"[step 2/3] Using provided .fa: {out_fsa}", file=sys.stderr)
        else:
            print(f"[step 2/3] Copying provided .fa to output location ...", file=sys.stderr)
            shutil.copy2(fsa_source, out_fsa)
    else:
        fsa_preexisted = os.path.exists(out_fsa)
        if fsa_preexisted:
            print(f"[step 2/3] Skipping — {out_fsa} already exists.", file=sys.stderr)
        else:
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

    if args.preconvert_only:
        print("\nDone (preconvert only).", file=sys.stderr)
        print(f"  Feature table : {out_tbl}", file=sys.stderr)
        print(f"  FASTA         : {out_fsa}", file=sys.stderr)
        return

    # ── Step 3: DDBJ MSS annotation ──────────────────────────────────────
    print("[step 3/3] Converting to DDBJ MSS annotation ...", file=sys.stderr)
    common: Optional[CommonModel] = None
    if args.common:
        try:
            common = load_common_json(args.common)
        except Exception as exc:
            sys.exit(f"Error: failed to load --common file '{args.common}':\n{exc}")

    sequence_roles = None
    if args.sequence_roles:
        try:
            sequence_roles = load_sequence_roles(args.sequence_roles)
        except Exception as exc:
            sys.exit(f"Error: failed to load sequence role file '{args.sequence_roles}':\n{exc}")

    write_ddbj_ann(
        out_tbl, out_fsa, out_ann,
        common=common,
        sequence_roles=sequence_roles,
        submission_category=args.submission_category,
    )

    if not args.keep_tmp and not tbl_preexisted:
        os.unlink(out_tbl)

    print("\nDone.", file=sys.stderr)
    print(f"  Annotation : {out_ann}", file=sys.stderr)
    print(f"  FASTA      : {out_fsa}", file=sys.stderr)
