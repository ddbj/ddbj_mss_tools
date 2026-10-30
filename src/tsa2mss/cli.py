"""
tsa2mss — Create DDBJ MSS TSA (Transcriptome Shotgun Assembly) submission files
from a transcript assembly FASTA.

Only TSA citing READ ARCHIVE (DBLINK "sequence read archive") is supported;
EST-based TSA with PRIMARY lines is not.

- DIVISION TSA, KEYWORD "TSA" / "Transcriptome Shotgun Assembly" are added
- source (in COMMON, 1..E) defaults to mol_type "transcribed RNA"
- required: DBLINK sequence read archive, SOURCE geo_loc_name / collection_date,
  ST_COMMENT (Assembly-Data) Assembly Method / Sequencing Technology
- --meta (metatranscriptome): KEYWORD ENV and /environmental_sample are added;
  SOURCE isolation_source is required. DIVISION stays TSA.

Missing required values are written empty with a warning.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from common.cli_args import add_output_args, validate_prefix, resolve_output
from common.fasta import write_clean_fasta
from common.gap_annotator import GapAnnotator
from common.models import load_common_json
from mss_builder.ann_writer import write_mss_ann


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="tsa2mss",
        description="Create DDBJ MSS TSA submission files from a transcript assembly FASTA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input", help="Input FASTA file (.fa / .fasta)")
    add_output_args(parser)
    parser.add_argument(
        "--common", required=True,
        help="JSON file with common submission metadata (DBLINK, SUBMITTER, REFERENCE, ST_COMMENT, SOURCE, …)",
    )
    parser.add_argument(
        "--meta", action="store_true",
        help="Metatranscriptome assembly (KEYWORD ENV, /environmental_sample, isolation_source required)",
    )
    args = parser.parse_args()

    validate_prefix(args, parser)
    in_path = Path(args.input).resolve()
    out_prefix, out_dir = resolve_output(args, in_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fsa = str(out_dir / out_prefix) + ".fa"
    out_ann = str(out_dir / out_prefix) + ".ann"

    try:
        common = load_common_json(args.common)
    except Exception as exc:
        sys.exit(f"Error: failed to load --common file '{args.common}':\n{exc}")

    common.SOURCE = {"mol_type": "transcribed RNA", **(common.SOURCE or {})}

    # TSA (non-CON) entries cannot carry assembly_gap /estimated_length="unknown"
    gaps = common.ASSEMBLY_GAP or []
    for cfg in gaps if isinstance(gaps, list) else [gaps]:
        if not cfg.enabled:
            continue
        resolved = GapAnnotator(
            linkage_evidence=cfg.linkage_evidence,
            gap_type=cfg.gap_type,
            estimated_length=cfg.estimated_length,
        ).estimated_length
        if resolved == "unknown":
            sys.exit(
                f"Error: ASSEMBLY_GAP (linkage_evidence={cfg.linkage_evidence!r}) resolves to "
                "estimated_length 'unknown', which TSA does not allow. Use 'paired-ends' or "
                "set estimated_length to 'known'."
            )

    print("[step 1/2] Writing clean FASTA with '//' separators ...", file=sys.stderr)
    if str(in_path) == out_fsa:
        sys.exit("Error: output FASTA would overwrite the input; use -o or -p.")
    write_clean_fasta(str(in_path), out_fsa)

    print("[step 2/2] Building DDBJ MSS annotation ...", file=sys.stderr)
    write_mss_ann(out_fsa, out_ann, common=common,
                  submission_category="TSA-ENV" if args.meta else "TSA")

    print("\nDone.", file=sys.stderr)
    print(f"  Annotation : {out_ann}", file=sys.stderr)
    print(f"  FASTA      : {out_fsa}", file=sys.stderr)
