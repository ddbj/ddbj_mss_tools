"""FASTA utilities shared by egapx2mss and wgs_maker."""

from __future__ import annotations

import gzip


def parse_fasta_lengths(fasta_path: str) -> dict[str, int]:
    """Return {seq_id: length} for every record in *fasta_path*."""
    return {seq_id: len(seq) for seq_id, seq in parse_fasta_sequences(fasta_path).items()}


def parse_fasta_sequences(fasta_path: str) -> dict[str, str]:
    """Return {seq_id: sequence} for every record in *fasta_path*."""
    with open(fasta_path) as fh:
        content = fh.read()

    seqs: dict[str, str] = {}
    blocks = content.split("\n>")
    for i, block in enumerate(blocks):
        if i == 0:
            if block.startswith(">"):
                block = block[1:]
            else:
                continue
        lines = block.splitlines()
        if not lines:
            continue
        header = lines[0].split()[0]
        if header.startswith("lcl|"):
            header = header[4:]
        seq = "".join(lines[1:]).rstrip("/")
        if header:
            seqs[header] = seq
    return seqs


def write_clean_fasta(raw_fsa: str, out_fsa: str) -> None:
    """
    Write FASTA with cleaned headers and DDBJ MSS separators.

    - Strip the 'lcl|' prefix from sequence IDs
    - Keep only the first whitespace-delimited token as the ID
    - Strip any trailing '/' or '//' separator from sequence content
    - Append '//' after each sequence entry on its own line
    """
    with open(out_fsa, "w") as fout:
        for seq_id, seq in parse_fasta_sequences(raw_fsa).items():
            seq_out = "\n".join(seq[i:i+70] for i in range(0, len(seq), 70))
            fout.write(f">{seq_id}\n{seq_out}\n//\n" if seq_out else f">{seq_id}\n//\n")


def read_fasta(file_name: str) -> list:
    """Read a FASTA file (optionally gzip-compressed) and return a list of SeqRecord objects."""
    from Bio import SeqIO
    if file_name.endswith(".gz"):
        return list(SeqIO.parse(gzip.open(file_name, "rt"), "fasta"))
    else:
        return list(SeqIO.parse(file_name, "fasta"))
