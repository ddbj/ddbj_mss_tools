from __future__ import annotations
import argparse
from gff2mss.assemble import build_ann_text
from gff2mss.emit import emit_fasta


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="gff2mss",
                                 description="Convert canonical INSDC GFF3 + FASTA to DDBJ MSS (.ann/.fasta)")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--mss-config", required=True)
    ap.add_argument("--common", required=True)
    ap.add_argument("--sequence-roles")
    ap.add_argument("--submission-category", default="")
    ap.add_argument("--locus-tag-start", type=int, default=None,
                    help="override [locus_tag].start (e.g. continue organelle numbering after nuclear)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    ann_text, seqs = build_ann_text(args.gff, args.fasta, args.mss_config, args.common,
                                    args.sequence_roles, args.submission_category,
                                    locus_tag_start=args.locus_tag_start)
    with open(f"{args.out}.ann", "w", encoding="utf-8") as fh:
        fh.write(ann_text)
    with open(f"{args.out}.fasta", "w", encoding="utf-8") as fh:
        fh.write(emit_fasta(seqs))
    print(f"[gff2mss] -> {args.out}.ann / {args.out}.fasta")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
