"""FASTA utilities shared by egapx2mss and wgs_maker."""

from __future__ import annotations

import gzip


def parse_fasta_lengths(fasta_path: str) -> dict[str, int]:
    """Return {seq_id: length} for every record in *fasta_path*."""
    lengths: dict[str, int] = {}
    current_id: str | None = None
    current_len = 0

    with open(fasta_path) as fh:
        for raw in fh:
            line = raw.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                header = line[1:].split()[0]
                if header.startswith("lcl|"):
                    header = header[4:]
                current_id = header
                current_len = 0
            elif current_id is not None:
                current_len += len(line)

    if current_id is not None:
        lengths[current_id] = current_len

    return lengths


def parse_fasta_sequences(fasta_path: str) -> dict[str, str]:
    """Return {seq_id: sequence} for every record in *fasta_path*."""
    seqs: dict[str, str] = {}
    current_id: str | None = None
    buf: list[str] = []

    with open(fasta_path) as fh:
        for raw in fh:
            line = raw.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(buf)
                header = line[1:].split()[0]
                if header.startswith("lcl|"):
                    header = header[4:]
                current_id = header
                buf = []
            elif current_id is not None:
                buf.append(line)

    if current_id is not None:
        seqs[current_id] = "".join(buf)

    return seqs


def write_clean_fasta(raw_fsa: str, out_fsa: str) -> None:
    """
    Write FASTA with cleaned headers and DDBJ MSS separators.

    - Strip the 'lcl|' prefix from sequence IDs
    - Keep only the first whitespace-delimited token as the ID
    - Append '//' after each sequence entry (DDBJ MSS requires this separator)

    Each sequence line is 70 bases and already ends with '\\n', so the
    separator is written as '//\\n' with no extra blank line in between,
    regardless of whether the total sequence length is a multiple of 70.
    """
    with open(raw_fsa) as fin, open(out_fsa, "w") as fout:
        in_seq = False
        for line in fin:
            if line.startswith(">"):
                if in_seq:
                    fout.write("//\n")
                seq_id = line[1:].split()[0]
                if seq_id.startswith("lcl|"):
                    seq_id = seq_id[4:]
                fout.write(f">{seq_id}\n")
                in_seq = True
            else:
                fout.write(line)
        if in_seq:
            fout.write("//\n")


def read_fasta(file_name: str) -> list:
    """Read a FASTA file (optionally gzip-compressed) and return a list of SeqRecord objects."""
    from Bio import SeqIO
    if file_name.endswith(".gz"):
        return list(SeqIO.parse(gzip.open(file_name, "rt"), "fasta"))
    else:
        return list(SeqIO.parse(file_name, "fasta"))
