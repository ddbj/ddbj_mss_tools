import os, tempfile
from Bio import SeqIO
from ddbj_gff import parse
from ddbj_gff.flatfile import flatfile_to_gff
from ddbj_gff.writer import write
from gff2mss.convert import convert
from gff2mss.emit import emit_ann
from gff2mss.config import MssConfig
import mss2ff.cli

FIX = os.path.join(os.path.dirname(__file__), "flatfile_fixtures", "trans_splicing_rps12.gbk")
PROT = ("MPTIQQLIRNKRQPIENRTKSPALKGCPQRRGVCTRVYTTTPKKPNSALRKIARVRLTSGF"
        "EITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIIRGTLDAVGVKDRQQGRSKYGVKKSK")


def test_trans_splicing_flatfile_roundtrip(tmp_path):
    rec = SeqIO.read(FIX, "genbank")
    # forward: flatfile -> canonical GFF -> gff2mss (.ann), organelle emit_mrna=false
    doc = parse(write(flatfile_to_gff(rec)))
    cfg = MssConfig(source={}, transl_table=11, product_default="hypothetical protein")
    cfg.emit_mrna = False
    mss_doc, _ = convert(doc, {rec.id: rec.seq}, cfg, common_rows=[])
    ann = emit_ann(mss_doc)

    # the .ann carries the trans CDS (join complement) + /trans_splicing + both introns
    assert "join(complement(1641..1754),93..324,829..854)" in ann
    assert "join(complement(855..1640),1..92)" in ann
    assert "325..828" in ann                     # cis intron
    assert "\ttrans_splicing\t" in ann
    assert ann.count("\tintron\t") == 2

    # reverse: .ann + fasta -> mss2ff -> flatfile'; CDS translation matches original
    ann_p = tmp_path / "o.ann"; ann_p.write_text(ann, encoding="utf-8")
    fa_p = tmp_path / "o.fasta"
    fa_p.write_text(f">{rec.id}\n{str(rec.seq)}\n", encoding="utf-8")
    ff_p = tmp_path / "o.ff"
    mss2ff.cli.main([str(ann_p), str(fa_p), "-o", str(ff_p)])
    rt = SeqIO.read(str(ff_p), "genbank")
    cds = [f for f in rt.features if f.type == "CDS"]
    assert len(cds) == 1
    prot = str(cds[0].extract(rt.seq).translate(table=11))
    prot = prot[:-1] if prot.endswith("*") else prot
    assert prot == PROT
    assert "trans_splicing" in cds[0].qualifiers
    assert sum(1 for f in rt.features if f.type == "intron") == 2
