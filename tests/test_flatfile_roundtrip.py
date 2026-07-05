import os
from Bio import SeqIO
from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.flatfile import flatfile_to_gff
from ddbj_gff.writer import write
from gff2mss.convert import convert
from gff2mss.emit import emit_ann
from gff2mss.config import MssConfig

FIX = os.path.join(os.path.dirname(__file__), "flatfile_fixtures", "citrus_unshiu_excerpt.gbk")


def _translations_from_flatfile(rec):
    out = []
    for f in rec.features:
        if f.type == "CDS":
            tt = int(f.qualifiers.get("transl_table", ["1"])[0])
            prot = str(f.extract(rec.seq).translate(table=tt))
            out.append(prot[:-1] if prot.endswith("*") else prot)
    return sorted(out)


def test_flatfile_to_gff_roundtrip_cds_translation():
    rec = SeqIO.read(FIX, "genbank")
    # forward: flatfile -> canonical GFF
    gff_text = write(flatfile_to_gff(rec))
    doc = parse(gff_text)
    seqs = {rec.id: rec.seq}
    cfg = MssConfig(source={}, transl_table=1, product_default="hypothetical protein")
    cfg.emit_mrna = True                                   # nuclear 3-level
    mss_doc, _ = convert(doc, seqs, cfg, common_rows=[])
    ann = emit_ann(mss_doc)

    # the .ann (which mss2ff renders verbatim to a flatfile) carries mRNA + CDS
    assert "\tCDS\t" in ann and "\tmRNA\t" in ann

    # CDS translations survive the loop: recompute from the round-tripped CDS locations
    # by re-parsing the GFF CDS spans and translating against the original sequence.
    from ddbj_gff.model import Feature
    cds_feats = [f for f in doc.features if f.type == "CDS"]
    got = []
    for c in cds_feats:
        loc = c.to_biopython_location()
        prot = str(loc.extract(rec.seq).translate(table=1))
        got.append(prot[:-1] if prot.endswith("*") else prot)
    assert sorted(got) == _translations_from_flatfile(rec)   # 3 CDS, translations match
