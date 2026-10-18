from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.model import Feature, Span
from gff2mss.config import MssConfig
from gff2mss.convert import build_entry_features, build_gene_features
from gff2mss.locus_tag import LocusTagAssigner


def cfg():
    return MssConfig(source={}, transl_table=1, product_default="hypothetical protein", locus_tag_prefix="PFX")


def test_gff_level_misc_feature_child_is_emitted():
    gene = Feature("g1", "S", "gene", [Span("c1", 1, 9, "+")], {}, [])
    mrna = Feature("g1.t1", "S", "mRNA", [Span("c1", 1, 9, "+")], {}, [])
    misc = Feature("g1.t1.cds", "S", "misc_feature", [Span("c1", 1, 3, "+", 0), Span("c1", 7, 9, "+", 0)],
                   {"note": ["nonfunctional CDS: internal stop codon(s) detected in g1.t1.cds; not translated"]}, [])
    mrna.children = [misc]
    gene.children = [mrna]
    diags = []
    feats = build_gene_features(gene, "all", LocusTagAssigner.from_config(cfg()), Seq("ATGTAATAA"), cfg(), diags)
    assert [f.key for f in feats] == ["misc_feature"]
    q = {x.key: x.value for x in feats[0].qualifiers}
    assert feats[0].location == "join(1..3,7..9)"
    assert q["locus_tag"] == "PFX_000010"
    assert q["note"].startswith("nonfunctional CDS") or "submitter_gene_id: g1" in q["note"]
    assert any(x.key == "note" and x.value.startswith("nonfunctional CDS") for x in feats[0].qualifiers)
    assert any(x.key == "note" and "submitter_transcript_id: g1.t1" in x.value for x in feats[0].qualifiers)
    assert not any(d.code == "no-exon" for d in diags)


def test_ignored_non_gene_roots_are_reported():
    gff = ("##gff-version 3\n"
           "c1\tS\tmRNA\t1\t9\t.\t+\t.\tID=m1\nc1\tS\tCDS\t1\t9\t.\t+\t0\tID=m1.c;Parent=m1\n"
           "c1\tS\tmRNA\t11\t19\t.\t+\t.\tID=m2\nc1\tS\tCDS\t11\t19\t.\t+\t0\tID=m2.c;Parent=m2\n")
    diags = []
    per = build_entry_features(parse(gff), {"c1": Seq("ATGAAATAAAATGAAATAA")}, cfg(), diags)
    assert per["c1"] == []
    d = next(d for d in diags if d.code == "root-ignored")
    assert "2" in d.message and "mRNA" in d.message and d.severity.value == "WARNING"
