from ann_meta.model import Entry
from ann_meta.validate import Issue, validate

COMMON = {
    "DBLINK": {"project": "PRJDB1", "sample": "SAMD1"},
    "SUBMITTER": {"ab_name": ["Mishima,T."], "contact": "Taro Mishima"},
    "SOURCE": {"organism": "Escherichia coli", "mol_type": "genomic DNA"},
}


def _codes(issues: list[Issue]) -> list[str]:
    return [i.message for i in issues]


def _v(common=None, entries=None, *, category="GNM", fasta_ids=None, allowed_unknown=frozenset()):
    return validate(COMMON if common is None else common,
                    [Entry("s1", 100)] if entries is None else entries,
                    category=category, fasta_ids=fasta_ids, allowed_unknown=set(allowed_unknown))


def test_valid_input_has_no_issues():
    entries = [Entry("chr", 100, type="chromosome", seq_name="1", status="complete", topology="circular"),
               Entry("p1", 50, type="plasmid", seq_name="pA", status="complete", topology="circular",
                     overrides={"strain": "K-12"}, extra=[("note", "a"), ("note", "b")])]
    assert _v(entries=entries, fasta_ids=["chr", "p1"]) == []


def test_rejects_tab_in_value():
    common = {**COMMON, "SUBMITTER": {"ab_name": ["Mishima,T."], "contact": "Taro\tMishima"}}
    issues = _v(common=common)
    assert Issue("SUBMITTER.contact", None, "forbidden_char") in issues
    issues = _v(entries=[Entry("s1", 100, overrides={"strain": "a\nb"})])
    assert Issue("strain", "s1", "forbidden_char") in issues
    issues = _v(entries=[Entry("s\x00", 100)])
    assert "forbidden_char" in _codes(issues)
    issues = _v(common={**COMMON, "_extra": {"SUBMITTER": [["fax", "1\x7f"]]}})
    assert "forbidden_char" in _codes(issues)
    issues = _v(entries=[Entry("s1", 100, extra=[("my\tname", "x")])], allowed_unknown={"my\tname"})
    assert "forbidden_char" in _codes(issues)


def test_rejects_unknown_qualifier_unless_allowed():
    entries = [Entry("s1", 100, overrides={"sub_strain": "x"})]
    assert Issue("sub_strain", "s1", "unknown_qualifier") in _v(entries=entries)
    assert _v(entries=entries, allowed_unknown={"sub_strain"}) == []
    common = {**COMMON, "SOURCE": {**COMMON["SOURCE"], "sub_strain": "x"}}
    assert Issue("SOURCE.sub_strain", None, "unknown_qualifier") in _v(common=common)
    assert _v(common=common, allowed_unknown={"sub_strain"}) == []
    # extra names go through the same gate: they do not widen the whitelist
    entries = [Entry("s1", 100, extra=[("my_flag", True)])]
    assert _v(entries=entries) == [Issue("my_flag", "s1", "unknown_qualifier")]
    assert _v(entries=entries, allowed_unknown={"my_flag"}) == []
    # an unknown extra name in one entry does not let SOURCE use it
    entries = [Entry("s1", 100, extra=[("sub_strain", "y")])]
    common = {**COMMON, "SOURCE": {**COMMON["SOURCE"], "sub_strain": "x"}}
    assert set(_codes(_v(common=common, entries=entries))) == {"unknown_qualifier"}
    assert len(_v(common=common, entries=entries)) == 2
    # role names left in extra by extract are known qualifiers
    assert _v(entries=[Entry("s1", 100, type="chromosome", seq_name="1", status="complete",
                             extra=[("plasmid", "p1")])]) == []


def test_tsv_typo_column_is_rejected(tmp_path):
    from ann_meta.entries_tsv import read_entries
    p = tmp_path / "e.tsv"
    p.write_text("#seq_id\ttype\tseq_name\tstatus\ttopology\tff_definition\tstrian\n"
                 "s1\tunplaced\t\tpartial\t\t\tK-12\n", encoding="utf-8")
    entries = read_entries(p)
    assert _v(entries=entries) == [Issue("strian", "s1", "unknown_qualifier")]


def test_rejects_line_separators_and_c1_controls():
    for bad in ("Hanako\u2028s1\tDATE", "a\u2028b", "a\u2029b", "a\x85b", "a\x1cb", "a\x0bb", "a\x9fb"):
        common = {**COMMON, "SUBMITTER": {"ab_name": ["Mishima,T."], "contact": bad}}
        assert Issue("SUBMITTER.contact", None, "forbidden_char") in _v(common=common), repr(bad)
        issues = _v(entries=[Entry("s1", 100, overrides={"strain": bad})])
        assert Issue("strain", "s1", "forbidden_char") in issues, repr(bad)
        assert "forbidden_char" in _codes(_v(entries=[Entry("s1", 100, ff_definition=bad)])), repr(bad)
    ok = {**COMMON, "SUBMITTER": {"ab_name": ["Müller,J."], "contact": "谷澤 靖洋 — Ångström café"}}
    assert _v(common=ok) == []
    assert _v(entries=[Entry("s1", 100, overrides={"strain": "株 NIG-1 µ"})]) == []


def test_fasta_id_missing_from_table():
    issues = _v(entries=[Entry("s1", 100)], fasta_ids=["s1", "s2"])
    assert issues == [Issue("entries", "s2", "seq_id_missing_from_table")]


def test_reserved_qualifier_in_overrides():
    issues = _v(entries=[Entry("s1", 100, overrides={"submitter_seqid": "x"})])
    assert Issue("submitter_seqid", "s1", "reserved_qualifier") in issues


def test_type_status_topology_domain():
    issues = _v(entries=[Entry("s1", 100, type="contig", status="done", topology="ring")])
    assert {"bad_type", "bad_status", "bad_topology"} <= set(_codes(issues))


def test_unplaced_and_partial_rules():
    issues = _v(entries=[Entry("s1", 100, seq_name="x", status="complete", topology="circular")])
    assert sorted(_codes(issues)) == ["unplaced_circular", "unplaced_not_partial", "unplaced_with_name"]
    issues = _v(entries=[Entry("c", 100, type="chromosome", seq_name="1", status="partial",
                               topology="circular")])
    assert Issue("topology", "c", "partial_circular") in issues


def test_plasmid_needs_name():
    issues = _v(entries=[Entry("p", 100, type="plasmid", status="complete")])
    assert Issue("seq_name", "p", "plasmid_needs_name") in issues
    issues = _v(entries=[Entry("o", 100, type="organelle", status="complete")])
    assert Issue("seq_name", "o", "organelle_needs_name") in issues


def test_two_chromosomes_need_names():
    one = [Entry("c1", 100, type="chromosome", status="complete")]
    assert _v(entries=one) == []
    two = [Entry("c1", 100, type="chromosome", seq_name="1", status="complete"),
           Entry("c2", 100, type="chromosome", status="complete")]
    assert _v(entries=two) == [Issue("seq_name", "c2", "chromosome_needs_name")]
    segs = [Entry("g1", 100, type="segment", status="complete"),
            Entry("g2", 100, type="segment", status="complete")]
    assert _codes(_v(entries=segs)) == ["segment_needs_name", "segment_needs_name"]


def test_role_not_allowed_for_wgs():
    entries = [Entry("c1", 100, type="chromosome", seq_name="1", status="complete"),
               Entry("s2", 100)]
    issues = _v(entries=entries, category="WGS")
    assert issues == [Issue("type", "c1", "role_not_allowed_for_category")]
    assert "role_not_allowed_for_category" in _codes(_v(entries=entries, category="GENE"))
    assert "role_not_allowed_for_category" in _codes(_v(entries=entries, category="TSA"))
    assert _v(entries=entries, category="GNM") == []
    assert _v(entries=entries, category="MAG") == []
    assert _v(entries=entries, category="NONE") == []


def test_seq_id_not_in_fasta():
    issues = _v(entries=[Entry("s1", 100), Entry("s2", 100)], fasta_ids=["s1"])
    assert issues == [Issue("seq_id", "s2", "seq_id_not_in_fasta")]
    assert _v(entries=[Entry("s1", 100), Entry("s2", 100)], fasta_ids=None) == []


def test_duplicate_seq_id():
    issues = _v(entries=[Entry("s1", 100), Entry("s1", 100)])
    assert issues == [Issue("seq_id", "s1", "duplicate_seq_id")]


def test_empty_seq_id_and_no_entries():
    assert Issue("seq_id", None, "empty_seq_id") in _v(entries=[Entry("", 100)])
    assert _v(entries=[]) == [Issue("entries", None, "no_entries")]


def test_unknown_category():
    assert _v(category="FOO") == [Issue("category", None, "unknown_category")]
    for ok in ("GENE", "NONE", "GNM", "WGS", "ENV", "MAG", "MAG-WGS", "TSA", "TSA-ENV"):
        assert _v(category=ok) == [], ok
