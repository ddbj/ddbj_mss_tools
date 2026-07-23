from __future__ import annotations
import argparse
import sys


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

    # gff2mss depends on the optional 'ddbj-gff' package. Import it lazily (after
    # argument parsing, so --help works without it) and fail with a clear message
    # when it is not installed — the other tools in this suite do not need it.
    try:
        from gff2mss.assemble import build_ann_text
        from gff2mss.emit import emit_fasta
    except ImportError as exc:
        print(
            "gff2mss requires the optional 'ddbj-gff' dependency, which is not installed.\n"
            "Install it with:  pip install 'ddbj-mss-tools[gff2mss]'\n"
            f"(import error: {exc})",
            file=sys.stderr,
        )
        return 1

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
