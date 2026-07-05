import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ddbj_gff import parse
from ddbj_gff.flatfile import flatfile_to_gff
from ddbj_gff.writer import write
from gff2mss.convert import convert
from gff2mss.emit import emit_ann
from gff2mss.config import MssConfig
import mss2ff.cli

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


def test_flatfile_to_gff_to_mss_to_flatfile_full_roundtrip(tmp_path):
    """Spec success-criterion (1): flatfile -> canonical GFF -> MSS .ann/.fsa
    -> mss2ff -> regenerated flatfile, and the regenerated CDS /translation
    values must match the original flatfile's CDS translations exactly.

    Unlike test_flatfile_to_gff_roundtrip_cds_translation above (which stops
    at the .ann and re-derives translations from GFF coordinates), this test
    actually runs mss2ff and inspects its output flatfile.
    """
    rec = SeqIO.read(FIX, "genbank")
    original_translations = _translations_from_flatfile(rec)

    # forward: flatfile -> canonical GFF -> MSS .ann
    gff_text = write(flatfile_to_gff(rec))
    doc = parse(gff_text)
    seqs = {rec.id: rec.seq}
    cfg = MssConfig(source={}, transl_table=1, product_default="hypothetical protein",
                     emit_mrna=True)
    mss_doc, _ = convert(doc, seqs, cfg, common_rows=[])
    ann_text = emit_ann(mss_doc)
    assert "\tCDS\t" in ann_text and "\tmRNA\t" in ann_text

    ann_path = tmp_path / "entry.ann"
    fasta_path = tmp_path / "entry.fasta"
    out_path = tmp_path / "entry.regenerated.gbk"
    ann_path.write_text(ann_text)
    SeqIO.write(SeqRecord(rec.seq, id=rec.id, description=""), str(fasta_path), "fasta")

    # reverse: mss2ff .ann + .fasta -> regenerated DDBJ flatfile
    mss2ff.cli.main([str(ann_path), str(fasta_path), "-o", str(out_path)])

    regenerated = SeqIO.read(str(out_path), "genbank")
    regenerated_cds = [f for f in regenerated.features if f.type == "CDS"]
    assert len(regenerated_cds) == 3

    # mss2ff computes and writes /translation directly onto each CDS
    # (see ff_writer._translate_cds) -- compare those, not a re-derivation.
    got = []
    for f in regenerated_cds:
        trans = f.qualifiers.get("translation")
        assert trans, f"regenerated CDS at {f.location} missing /translation"
        prot = trans[0]
        got.append(prot[:-1] if prot.endswith("*") else prot)

    assert sorted(got) == original_translations
