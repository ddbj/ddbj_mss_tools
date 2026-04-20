"""Command-line interface for mss2ff: MSS annotation → DDBJ flat file."""

from __future__ import annotations

import argparse
import sys
from datetime import date, datetime
from pathlib import Path

from Bio import SeqIO


def _parse_date(s: str) -> date:
    for fmt in ("%Y-%m-%d", "%d-%b-%Y", "%Y%m%d"):
        try:
            return datetime.strptime(s, fmt).date()
        except ValueError:
            continue
    raise argparse.ArgumentTypeError(
        f"Cannot parse date {s!r}. Use YYYY-MM-DD or DD-MON-YYYY."
    )


def _read_fasta(path: str) -> dict[str, str]:
    """Read FASTA/FSA file, returning {seq_id: sequence_str} (lowercase).

    Handles DDBJ FSA format where entries are separated by '//' lines.
    """
    seqs: dict[str, str] = {}
    cur_id: str | None = None
    cur_chunks: list[str] = []

    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_chunks).lower()
                cur_id = line[1:].split()[0]
                cur_chunks = []
            elif line.strip() == "//":
                continue  # skip FSA entry separators
            elif cur_id is not None:
                cur_chunks.append(line.strip())

    if cur_id is not None:
        seqs[cur_id] = "".join(cur_chunks).lower()

    return seqs


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="mss2ff",
        description="Convert DDBJ MSS annotation file to DDBJ flat file format.",
    )
    parser.add_argument(
        "ann",
        metavar="ANN",
        help="MSS annotation file (.ann or .annt.tsv)",
    )
    parser.add_argument(
        "--fasta", "-f",
        metavar="FASTA",
        help="FASTA sequence file (.fa / .fsa / .fasta). Required for translation.",
    )
    parser.add_argument(
        "--output", "-o",
        metavar="OUTPUT",
        default="-",
        help="Output flat file path (default: stdout)",
    )
    parser.add_argument(
        "--division", "-d",
        metavar="DIV",
        default="UNK",
        help="DDBJ division code (e.g. BCT, VRL, PLN). Default: UNK",
    )
    parser.add_argument(
        "--submission-date", "-s",
        metavar="DATE",
        type=_parse_date,
        default=None,
        help="Submission date for Reference 1 JOURNAL line (YYYY-MM-DD). Default: today",
    )
    parser.add_argument(
        "--file-date",
        metavar="DATE",
        type=_parse_date,
        default=None,
        help="File creation date for LOCUS line (YYYY-MM-DD). Default: today",
    )
    parser.add_argument(
        "--email",
        metavar="EMAIL",
        default="mss2ff@ddbj.nig.ac.jp",
        help="Email address for NCBI Entrez API calls",
    )
    parser.add_argument(
        "--no-taxonomy",
        action="store_true",
        help="Skip NCBI taxonomy lookup (no lineage or taxon db_xref)",
    )

    args = parser.parse_args(argv)

    # Lazy imports (keep startup fast)
    from .ann_parser import parse_ann
    from .ff_writer import write_ff

    # Parse annotation file
    try:
        common, entries = parse_ann(args.ann)
    except Exception as exc:
        print(f"mss2ff: Error reading annotation file: {exc}", file=sys.stderr)
        sys.exit(1)

    if not entries:
        print("mss2ff: No entries found in annotation file.", file=sys.stderr)
        sys.exit(1)

    # Read sequences
    sequences: dict[str, str] = {}
    if args.fasta:
        try:
            sequences = _read_fasta(args.fasta)
        except Exception as exc:
            print(f"mss2ff: Error reading FASTA file: {exc}", file=sys.stderr)
            sys.exit(1)

    # Open output
    if args.output == "-":
        out = sys.stdout
    else:
        try:
            out = open(args.output, "w", encoding="utf-8")
        except OSError as exc:
            print(f"mss2ff: Cannot open output file: {exc}", file=sys.stderr)
            sys.exit(1)

    try:
        write_ff(
            common=common,
            entries=entries,
            sequences=sequences,
            output=out,
            division=args.division,
            submission_date=args.submission_date,
            file_date=args.file_date,
            email=args.email,
            no_taxonomy=args.no_taxonomy,
        )
    finally:
        if args.output != "-":
            out.close()


if __name__ == "__main__":
    main()
