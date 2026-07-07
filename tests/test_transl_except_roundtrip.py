import os
from Bio import SeqIO
from ddbj_gff import parse
from ddbj_gff.flatfile import flatfile_to_gff
from ddbj_gff.writer import write
from gff2mss.convert import build_entry_features
from gff2mss.config import MssConfig

FIX = os.path.join(os.path.dirname(__file__), "flatfile_fixtures", "transl_except_p87.gbk")


def _cds_feats():
    rec = SeqIO.read(FIX, "genbank")
    doc = parse(write(flatfile_to_gff(rec)))               # flatfile -> canonical GFF (recoded_codon child)
    seqs = {rec.id: rec.seq}
    cfg = MssConfig(source={}, transl_table=1, product_default="hypothetical protein")
    cfg.emit_mrna = False                                   # virus 2-level
    per_entry = build_entry_features(doc, seqs, cfg, [])
    return [f for feats in per_entry.values() for f in feats]


def test_transl_except_qualifier_emitted_on_clean_cds():
    feats = _cds_feats()
    cds = [f for f in feats if f.key == "CDS"]
    assert len(cds) == 1                                    # translated cleanly via Pyl -> CDS, not misc_feature
    te = [q.value for q in cds[0].qualifiers if q.key == "transl_except"]
    assert te == ["(pos:746..748,aa:Pyl)"]                  # qualifier emitted, matches original
