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


def test_transl_except_full_roundtrip(tmp_path):
    from gff2mss.convert import convert
    from gff2mss.emit import emit_ann
    import mss2ff.cli
    rec = SeqIO.read(FIX, "genbank")
    doc = parse(write(flatfile_to_gff(rec)))
    cfg = MssConfig(source={}, transl_table=1, product_default="hypothetical protein")
    cfg.emit_mrna = False
    mss_doc, _ = convert(doc, {rec.id: rec.seq}, cfg, common_rows=[])
    ann = emit_ann(mss_doc)
    assert "\ttransl_except\t(pos:746..748,aa:Pyl)" in ann     # qualifier in the .ann

    ann_p = tmp_path / "o.ann"; ann_p.write_text(ann, encoding="utf-8")
    fa_p = tmp_path / "o.fasta"; fa_p.write_text(f">{rec.id}\n{str(rec.seq)}\n", encoding="utf-8")
    ff_p = tmp_path / "o.ff"
    mss2ff.cli.main([str(ann_p), str(fa_p), "-o", str(ff_p)])
    rt = SeqIO.read(str(ff_p), "genbank")
    cds = [f for f in rt.features if f.type == "CDS"]
    assert len(cds) == 1
    assert cds[0].qualifiers.get("transl_except") == ["(pos:746..748,aa:Pyl)"]   # /transl_except preserved
    assert "O" in cds[0].qualifiers["translation"][0]           # Pyl -> O in the regenerated protein
