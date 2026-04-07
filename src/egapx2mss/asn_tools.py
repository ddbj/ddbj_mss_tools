"""
Download and run NCBI command-line tools (asn2gb, asn2fsa).

Tools are downloaded from https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/cmdline/
and cached in a user-specified directory (default: ~/.local/share/ddbj_mss_tools/bin).

If a cached binary produces no valid output (indicating it has expired), it is
automatically re-downloaded and the operation is retried once.
"""

from __future__ import annotations

import gzip
import os
import platform
import shutil
import stat
import subprocess
import sys
import urllib.request
from pathlib import Path


NCBI_CMDLINE_URL = "https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/cmdline/"

_PLATFORM_SUFFIX: dict[str, dict[str, str]] = {
    "Darwin": {"asn2gb": "asn2gb.mac.gz",     "asn2fsa": "asn2fsa.mac.gz"},
    "Linux":  {"asn2gb": "asn2gb.linux64.gz", "asn2fsa": "asn2fsa.linux64.gz"},
}

DEFAULT_BIN_DIR = Path(__file__).parent.parent.parent / "bin"


# ── Download ──────────────────────────────────────────────────────────────────

def _download_tool(name: str, bin_dir: Path, force: bool = False) -> Path:
    """
    Download *name* (asn2gb or asn2fsa) to *bin_dir*.

    If *force* is False (default) and the binary already exists, return it as-is.
    If *force* is True, always fetch the latest version from NCBI (overwrites the
    existing file), which is used when the current binary has expired.
    """
    dest = bin_dir / name
    if dest.exists() and not force:
        return dest

    system = platform.system()
    if system not in _PLATFORM_SUFFIX:
        raise RuntimeError(
            f"Unsupported platform '{system}'. "
            f"Please download {name} manually from {NCBI_CMDLINE_URL}"
        )

    bin_dir.mkdir(parents=True, exist_ok=True)
    gz_name = _PLATFORM_SUFFIX[system][name]
    url = NCBI_CMDLINE_URL + gz_name
    gz_path = bin_dir / gz_name

    action = "Re-downloading" if (dest.exists() and force) else "Downloading"
    print(f"[setup] {action} {name} from {url} ...", file=sys.stderr)
    urllib.request.urlretrieve(url, gz_path)

    with gzip.open(gz_path, "rb") as f_in, open(dest, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    dest.chmod(dest.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    gz_path.unlink()
    print(f"[setup] Saved to {dest}", file=sys.stderr)
    return dest


def ensure_tools(bin_dir: Path) -> tuple[Path, Path]:
    """Return paths to asn2gb and asn2fsa, downloading them if necessary."""
    return _download_tool("asn2gb", bin_dir), _download_tool("asn2fsa", bin_dir)


# ── Output validation ─────────────────────────────────────────────────────────

def _is_valid_tbl(path: str) -> bool:
    """A valid feature table contains at least one '>Feature' header line."""
    try:
        with open(path) as fh:
            return any(line.startswith(">Feature") for line in fh)
    except OSError:
        return False


def _is_valid_fsa(path: str) -> bool:
    """A valid FASTA file contains at least one '>' header line."""
    try:
        with open(path) as fh:
            return any(line.startswith(">") for line in fh)
    except OSError:
        return False


# ── ASN.1 block splitting ────────────────────────────────────────────────────

def _iter_asn_blocks(filepath: str):
    """Yield line lists for each top-level Seq-entry block in a catenated ASN.1 file."""
    block: list[str] = []
    with open(filepath) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("Seq-entry ::=") and block:
                yield block
                block = []
            block.append(line)
    if block:
        yield block


# ── Tool runners ──────────────────────────────────────────────────────────────

def _run_asn2gb_tbl_once(asn2gb: Path, asn_path: str, out_tbl: str) -> None:
    with open(out_tbl, "w") as fout:
        subprocess.run(
            [str(asn2gb), "-f", "t", "-a", "q", "-i", asn_path],
            stdout=fout,
            stderr=subprocess.PIPE,
            check=True,
        )


def run_asn2gb_tbl(asn2gb: Path, asn_path: str, out_tbl: str) -> Path:
    """
    Run asn2gb in feature-table mode (-f t) on a catenated ASN.1 file.

    If the output is empty or contains no '>Feature' headers (which happens when
    the binary has expired), the tool is re-downloaded and the run is retried
    once before raising an error.

    Returns the (possibly updated) path to the asn2gb binary.
    """
    _run_asn2gb_tbl_once(asn2gb, asn_path, out_tbl)

    if not _is_valid_tbl(out_tbl):
        print(
            "[asn2gb]  Output is empty or invalid — binary may have expired. "
            "Re-downloading ...",
            file=sys.stderr,
        )
        asn2gb = _download_tool("asn2gb", asn2gb.parent, force=True)
        _run_asn2gb_tbl_once(asn2gb, asn_path, out_tbl)
        if not _is_valid_tbl(out_tbl):
            raise RuntimeError(
                "asn2gb produced no valid output even after re-download. "
                "Check that the input file is a valid ASN.1 file."
            )

    print(f"[asn2gb]  → {out_tbl}", file=sys.stderr)
    return asn2gb


def _run_asn2fsa_once(asn2fsa: Path, asn_path: str, tmpdir: str) -> list[str]:
    parts: list[str] = []

    for i, lines in enumerate(_iter_asn_blocks(asn_path)):
        tmp_asn = os.path.join(tmpdir, f"record_{i}.asn")
        tmp_fsa = os.path.join(tmpdir, f"record_{i}.fsa")

        with open(tmp_asn, "w") as f:
            f.write("\n".join(lines) + "\n")

        subprocess.run(
            [str(asn2fsa), "-a", "a", "-i", tmp_asn, "-o", tmp_fsa],
            capture_output=True,
        )

        if os.path.exists(tmp_fsa):
            content = open(tmp_fsa).read()
            if content.strip():
                parts.append(content)

    return parts


def run_asn2fsa(asn2fsa: Path, asn_path: str, out_fsa: str, tmpdir: str) -> Path:
    """
    Run asn2fsa on every Seq-entry block in a catenated ASN.1 file.

    asn2fsa does not support catenated Seq-entry files natively, so each block
    is written to a temporary file and processed individually; the results are
    concatenated into *out_fsa*.

    If the combined output contains no FASTA sequences (which happens when the
    binary has expired), the tool is re-downloaded and the run is retried once.

    Returns the (possibly updated) path to the asn2fsa binary.
    """
    parts = _run_asn2fsa_once(asn2fsa, asn_path, tmpdir)

    if not parts:
        print(
            "[asn2fsa] Output is empty — binary may have expired. "
            "Re-downloading ...",
            file=sys.stderr,
        )
        asn2fsa = _download_tool("asn2fsa", asn2fsa.parent, force=True)
        parts = _run_asn2fsa_once(asn2fsa, asn_path, tmpdir)
        if not parts:
            raise RuntimeError(
                "asn2fsa produced no output even after re-download. "
                "Check that the input file is a valid ASN.1 file."
            )

    with open(out_fsa, "w") as f:
        f.write("".join(parts))

    print(f"[asn2fsa] → {out_fsa}  ({len(parts)} sequences)", file=sys.stderr)
    return asn2fsa
