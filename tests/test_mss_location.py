from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from gff2mss.config import MssConfig
from gff2mss.convert import collect_spans, build_insdc_location, extract_seq, build_cds_feature


def test_collect_spans_unions_children_of_type():
    mrna = Feature("m", "S", "mRNA", [Span("c", 1, 100, "+")], {}, [])
    e1 = Feature("e1", "S", "exon", [Span("c", 1, 10, "+")], {}, [])
    e2 = Feature("e2", "S", "exon", [Span("c", 20, 30, "+")], {}, [])
    c1 = Feature("cds", "S", "CDS", [Span("c", 3, 10, "+"), Span("c", 20, 28, "+")], {}, [])
    mrna.children = [e1, e2, c1]
    assert len(collect_spans(mrna, "exon")) == 2
    assert len(collect_spans(mrna, "CDS")) == 2  # NCBI-style single feature with 2 spans


def test_plus_strand_join_and_partials():
    spans = [Span("c", 20, 30, "+"), Span("c", 1, 10, "+")]  # unsorted input
    assert build_insdc_location(spans, 10000) == "join(1..10,20..30)"
    assert build_insdc_location([Span("c", 1, 10, "+")], 10000, five_prime_partial=True) == "<1..10"
    assert build_insdc_location([Span("c", 1, 10, "+")], 10000, three_prime_partial=True) == "1..>10"
    assert build_insdc_location(spans, 10000, five_prime_partial=True, three_prime_partial=True) == "join(<1..10,20..>30)"


def test_minus_strand_complement_ascending_inside():
    spans = [Span("c", 1, 10, "-"), Span("c", 20, 30, "-")]
    assert build_insdc_location(spans, 10000) == "complement(join(1..10,20..30))"
    assert build_insdc_location([Span("c", 5, 9, "-")], 10000) == "complement(5..9)"
    # 5' end of a minus feature is the high-coordinate end -> AfterPosition
    assert build_insdc_location([Span("c", 5, 9, "-")], 10000, five_prime_partial=True) == "complement(5..>9)"
    assert build_insdc_location([Span("c", 5, 9, "-")], 10000, three_prime_partial=True) == "complement(<5..9)"


def test_extract_plus_and_minus():
    genome = Seq("ATGAAATAA")            # +: ATGAAATAA
    assert str(extract_seq([Span("c", 1, 9, "+")], genome)) == "ATGAAATAA"
    rc = Seq("TTATTTCAT")                # revcomp of ATGAAATAA
    assert str(extract_seq([Span("c", 1, 9, "-")], rc)) == "ATGAAATAA"


def test_minus_strand_origin_spanning():
    assert build_insdc_location([Span("c", 4447, 5268, "-")], 5125) == "complement(join(4447..5125,1..143))"


def test_plus_strand_origin_spanning():
    assert build_insdc_location([Span("c", 4447, 5268, "+")], 5125) == "join(4447..5125,1..143)"


def test_origin_spanning_partials_minus():
    assert build_insdc_location([Span("c", 4447, 5268, "-")], 5125,
                                five_prime_partial=True) == "complement(join(4447..5125,1..>143))"
    assert build_insdc_location([Span("c", 4447, 5268, "-")], 5125,
                                three_prime_partial=True) == "complement(join(<4447..5125,1..143))"


def test_extract_origin_spanning_plus_translates():
    # 9 bp circular; plus CDS 7..12 wraps: head=7..9 "ATG", tail=1..3 "TAA" -> "ATGTAA" -> M*
    genome = Seq("TAACCCATG")
    ex = extract_seq([Span("c", 7, 12, "+")], genome)
    assert str(ex) == "ATGTAA"
    assert str(ex.translate(table=11)) == "M*"


def test_extract_origin_spanning_minus_translates():
    # minus CDS 7..12 on revcomp genome yields the same coding sequence
    genome = Seq("TAACCCATG").reverse_complement()  # so minus strand of 7..12 -> ATGTAA
    ex = extract_seq([Span("c", 7, 12, "-")], genome)
    assert str(ex.translate(table=11)) == "M*"


def test_multi_exon_wrap_warns():
    # NOTE: the task brief's sketch called this with an undefined `_MinimalCfg()`; no such
    # helper exists in this repo. Following the brief's own contingency ("confirm the exact
    # build_cds_feature signature ... do not invent parameters"), this uses the same
    # MssConfig(...) construction already established by tests/test_mss_cds.py::cfg().
    cds = Feature("cds", "S", "CDS",
                  [Span("c", 1, 6, "+"), Span("c", 7, 14, "+")], {"transl_table": ["11"]}, [])  # 14 > 10
    mrna = Feature("m", "S", "mRNA", [Span("c", 1, 14, "+")], {}, [])
    mrna.children = [cds]
    gene = Feature("g", "S", "gene", [Span("c", 1, 14, "+")], {}, [])
    cfg = MssConfig(source={}, transl_table=1, product_default="hypothetical protein")
    diags = []
    build_cds_feature(mrna, gene, "T_0001", Seq("ATGAAATAACC"[:10] + "A"), cfg, diags)
    assert any(d.code == "multi-exon-origin-spanning" for d in diags)


def test_transl_except_origin_spanning_warns():
    # 9 bp circular genome; CDS 7..12 wraps (head=7..9 "ATG", tail=1..3 "TAA") and
    # carries a raw transl_except. build_cds_feature's transl_except translation
    # path (SeqFeature built from cds_feat.to_biopython_location()) is NOT
    # wrap-aware, so this combination must raise a diagnostic warning the
    # submitter that the resulting protein may be wrong.
    genome = Seq("TAACCCATG")
    cds = Feature("cds", "S", "CDS", [Span("c", 7, 12, "+", 0)],
                  {"transl_except": ["(pos:7,aa:Term)"]}, [])
    mrna = Feature("m", "S", "mRNA", [Span("c", 7, 12, "+")], {}, [])
    mrna.children = [cds]
    gene = Feature("g", "S", "gene", [Span("c", 7, 12, "+")], {}, [])
    cfg = MssConfig(source={}, transl_table=1, product_default="hypothetical protein")
    diags = []
    build_cds_feature(mrna, gene, "T_0001", genome, cfg, diags)
    assert any(d.code == "transl-except-origin-spanning" for d in diags)
