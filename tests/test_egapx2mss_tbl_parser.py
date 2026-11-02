from egapx2mss.tbl_parser import collect_qualifiers, format_location


def test_minus_strand_partials():
    # NCBI tbl: minus strand lists the 5' end first; '<' = 5'-partial, '>' = 3'-partial
    assert format_location([("<300", "200"), ("150", "100")]) == "complement(join(100..150,200..>300))"
    assert format_location([("300", "200"), ("150", ">100")]) == "complement(join(<100..150,200..300))"
    assert format_location([("<100", "200"), ("300", ">400")]) == "join(<100..200,300..>400)"


def _quals(ftype, **q):
    return collect_qualifiers({"type": ftype, "intervals": [], "qualifiers": q})


def test_cds_keeps_transl_except():
    q = _quals("CDS", transl_except=["(pos:complement(10..12),aa:Sec)"])
    assert ("transl_except", "(pos:complement(10..12),aa:Sec)") in q


def test_low_quality_exception_becomes_artificial_location():
    for ftype in ("CDS", "mRNA"):
        q = _quals(ftype, exception=["low-quality sequence region"])
        assert ("artificial_location", "low-quality sequence region") in q
        assert not any(n == "exception" for n, _ in q)
    q = _quals("CDS", exception=["annotated by transcript or proteomic data"])
    assert ("exception", "annotated by transcript or proteomic data") in q


def test_standard_name_becomes_note():
    q = _quals("V_segment", locus_tag=["X_1"], standard_name=["T cell receptor alpha variable 9-2-like"])
    assert q == [("locus_tag", "X_1"), ("note", "T cell receptor alpha variable 9-2-like")]
