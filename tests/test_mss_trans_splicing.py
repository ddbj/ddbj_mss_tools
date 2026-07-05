import os
from Bio import SeqIO
from ddbj_gff import parse
from ddbj_gff.normalize.normalize import normalize
from gff2mss.convert import build_entry_features
from gff2mss.config import MssConfig

FIX = os.path.join(os.path.dirname(__file__), "mss_fixtures")
PROT = ("MPTIQQLIRNKRQPIENRTKSPALKGCPQRRGVCTRVYTTTPKKPNSALRKIARVRLTSGF"
        "EITAYIPGIGHNLQEHSVVLVRGGRVKDLPGVRYHIIRGTLDAVGVKDRQQGRSKYGVKKSK")


def _load():
    with open(os.path.join(FIX, "trans_splicing_rps12.gff3")) as fh:
        doc, _ = normalize(parse(fh.read()))
    seqs = {r.id: r.seq for r in SeqIO.parse(os.path.join(FIX, "trans_splicing_rps12.fasta"), "fasta")}
    cfg = MssConfig(source={}, transl_table=11, product_default="hypothetical protein")
    cfg.emit_mrna = False
    return doc, seqs, cfg


def test_trans_spliced_cds_location_and_translation():
    doc, seqs, cfg = _load()
    feats = build_entry_features(doc, seqs, cfg, [])["AP025455.1"]
    cds = [f for f in feats if f.key == "CDS"]
    assert len(cds) == 1
    assert cds[0].location == "join(complement(1641..1754),93..324,829..854)"
    # clean translation -> CDS (not misc_feature)
    assert not any(f.key == "misc_feature" for f in feats)
    # /trans_splicing qualifier present (valueless)
    assert any(q.key == "trans_splicing" for q in cds[0].qualifiers)
    # correct product + codon_start
    assert any(q.key == "product" and q.value == "ribosomal protein S12" for q in cds[0].qualifiers)
    assert any(q.key == "codon_start" and q.value == "1" for q in cds[0].qualifiers)


def test_intron_features_emitted():
    doc, seqs, cfg = _load()
    feats = build_entry_features(doc, seqs, cfg, [])["AP025455.1"]
    introns = [f for f in feats if f.key == "intron"]
    locs = {f.location for f in introns}
    assert "join(complement(855..1640),1..92)" in locs   # trans intron 1
    assert "325..828" in locs                             # cis intron 2
    trans_intron = next(f for f in introns if f.location == "join(complement(855..1640),1..92)")
    assert any(q.key == "trans_splicing" for q in trans_intron.qualifiers)
    assert any(q.key == "number" and q.value == "1" for q in trans_intron.qualifiers)
    assert any(q.key == "gene" and q.value == "rps12" for q in trans_intron.qualifiers)
    cis_intron = next(f for f in introns if f.location == "325..828")
    assert not any(q.key == "trans_splicing" for q in cis_intron.qualifiers)


def test_ann_text_end_to_end():
    from gff2mss.convert import convert
    from gff2mss.emit import emit_ann
    doc, seqs, cfg = _load()                     # parse + normalize (location= built)
    mss_doc, _ = convert(doc, seqs, cfg, common_rows=[])
    text = emit_ann(mss_doc)
    assert "join(complement(1641..1754),93..324,829..854)" in text   # trans CDS
    assert "join(complement(855..1640),1..92)" in text               # trans intron 1
    assert "\ttrans_splicing\t" in text          # valueless qualifier row (empty value col)
    assert "\tintron\t" in text                  # intron feature key column
