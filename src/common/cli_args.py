"""
Shared argparse helpers for output-path options used across MSS tools.

Provides a consistent -o/--outdir and -p/--prefix pattern:
  add_output_args()  — registers both arguments on a parser
  validate_prefix()  — validates that --prefix contains no path separators
  resolve_output()   — resolves (out_prefix, out_dir) from args + a reference path
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path


def add_output_args(parser: argparse.ArgumentParser) -> None:
    """Add -o/--outdir and -p/--prefix to *parser*."""
    parser.add_argument(
        "-o", "--outdir",
        default=None,
        help=(
            "Output directory (created if absent). "
            "Default: same directory as input file."
        ),
    )
    parser.add_argument(
        "-p", "--prefix",
        default=None,
        help=(
            "Output filename prefix (basename only, no directory separators). "
            "Default: basename of input file without extension."
        ),
    )


def validate_prefix(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
) -> None:
    """Call parser.error() if args.prefix contains a directory separator."""
    if args.prefix and (os.sep in args.prefix or "/" in args.prefix):
        parser.error(
            "--prefix must be a basename only (no directory separators); "
            "use --outdir for the directory"
        )


def resolve_output(
    args: argparse.Namespace,
    ref_path: Path,
) -> tuple[str, Path]:
    """
    Return (out_prefix, out_dir) derived from *args* and *ref_path*.

    out_prefix = args.prefix  if given, else ref_path.stem
    out_dir    = args.outdir  if given, else ref_path.parent

    The caller is responsible for calling out_dir.mkdir(parents=True, exist_ok=True).
    """
    out_prefix: str  = args.prefix if args.prefix else ref_path.stem
    out_dir: Path    = Path(args.outdir).resolve() if args.outdir else ref_path.parent
    return out_prefix, out_dir
