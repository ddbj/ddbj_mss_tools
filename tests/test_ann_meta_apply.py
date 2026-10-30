import dataclasses
from pathlib import Path

import pytest

from ann_fix.annfile import AnnFile
from ann_meta.apply import apply
from ann_meta.entries_tsv import read_entries, write_entries
from ann_meta.extract import COMMON_TARGETS, extract, read_fasta_lengths
from ann_meta.model import Entry, MetaDoc
from ann_meta.validate import Issue, validate

FIX = Path(__file__).parent / "ann_meta_fixtures"
FIXTURES = ("dfast_draft", "dfast_complete", "marker_set", "handmade")


# -- helpers -------------------------------------------------------------------

def _load(name: str) -> tuple[AnnFile, MetaDoc, dict[str, int]]:
    lengths = read_fasta_lengths(FIX / f"{name}.fasta")
    ann = AnnFile.read(FIX / f"{name}.ann")
    return ann, extract(ann, lengths), lengths


def _common(doc: MetaDoc) -> dict:
    """The common JSON as GGS assembles it from a MetaDoc."""
    c = dict(doc.common)
    c["SOURCE"] = dict(doc.source)
    if doc.source_identifier:
        c["SOURCE_IDENTIFIER"] = doc.source_identifier
    if doc.common_extra:
        c["_extra"] = {k: [list(t) for t in v] for k, v in doc.common_extra.items()}
    return c


def _category(doc: MetaDoc) -> str:
    return doc.category.value or "GENE"


def _apply_file(src: Path, out: Path, *, common=None, entries=None, category=None,
                lengths=None, mutate=None) -> Path:
    ann = AnnFile.read(src)
    doc = extract(ann, lengths)
    common = _common(doc) if common is None else common
    entries = doc.entries if entries is None else entries
    if mutate:
        mutate(common, entries)
    apply(ann, common, entries, category=category or _category(doc)).write(out)
    return out


def _feature_rows(path: Path) -> list[tuple[str, bytes, bool]]:
    """Every row outside source / TOPOLOGY / the 6 COMMON targets, in order.

    Each item is (entry, raw line incl. eol, first kept row of its entry).
    """
    out = []
    entry = key = ""
    seen: set[str] = set()
    for raw in path.read_bytes().splitlines(keepends=True):
        cols = raw.rstrip(b"\r\n").split(b"\t")
        cols += [b""] * (5 - len(cols))
        if cols[0]:
            entry = cols[0].decode()
        if cols[1]:
            key = cols[1].decode()
        if key in ("source", "TOPOLOGY") or (entry == "COMMON" and key in COMMON_TARGETS):
            continue
        out.append((entry, raw, entry not in seen))
        seen.add(entry)
    return out


def _assert_feature_rows_equal(src: Path, out: Path) -> None:
    """Kept rows are byte-identical (column 1 included).

    Only exceptions, both on the first remaining row of an entry, where the entry
    marker legitimately moves (columns 2-5 and the eol must still match):
    - the entry's source was removed (COMMON-source layout): the row takes over
      the marker, so column 1 must be the entry id;
    - a source was inserted for an entry that had none: the new source takes
      the marker, so column 1 must become empty.
    """
    a, b = _feature_rows(src), _feature_rows(out)
    assert [(e, first) for e, _, first in a] == [(e, first) for e, _, first in b]
    in_sources = {f.entry for f in AnnFile.read(src).features() if f.key == "source"}
    out_sources = {f.entry for f in AnnFile.read(out).features() if f.key == "source"}
    for (entry, raw_a, first), (_, raw_b, _) in zip(a, b):
        if raw_a == raw_b:
            continue
        assert first, (raw_a, raw_b)
        marker = entry.encode()
        if entry in in_sources and entry not in out_sources:
            assert raw_a.startswith(b"\t") and raw_b == marker + raw_a, (raw_a, raw_b)
        else:
            assert entry not in in_sources and entry in out_sources, (raw_a, raw_b)
            assert raw_b.startswith(b"\t") and raw_a == marker + raw_b, (raw_a, raw_b)


def _entry_markers(path: Path) -> list[str]:
    return [line.split(b"\t", 1)[0].decode() for line in path.read_bytes().splitlines()
            if not line.startswith(b"\t")]


def _no_bare_lf(data: bytes) -> bool:
    return data.count(b"\n") == data.count(b"\r\n")


# -- round trips ---------------------------------------------------------------

@pytest.mark.parametrize("name", FIXTURES)
def test_roundtrip_feature_rows_identical(name, tmp_path):
    src = FIX / f"{name}.ann"
    out = _apply_file(src, tmp_path / "out.ann", lengths=read_fasta_lengths(FIX / f"{name}.fasta"))
    _assert_feature_rows_equal(src, out)
    assert [m for m in _entry_markers(out)] == [m for m in _entry_markers(src)]
    if name == "handmade":
        assert _no_bare_lf(out.read_bytes())
    else:
        assert b"\r" not in out.read_bytes()


@pytest.mark.parametrize("name", FIXTURES)
def test_roundtrip_feature_rows_identical_crlf_and_no_final_newline(name, tmp_path):
    data = (FIX / f"{name}.ann").read_bytes().replace(b"\r\n", b"\n")
    for variant, body in (("crlf", data.replace(b"\n", b"\r\n")), ("lf_nofinal", data.rstrip(b"\n")),
                          ("crlf_nofinal", data.replace(b"\n", b"\r\n").rstrip(b"\r\n"))):
        src = tmp_path / f"{variant}.ann"
        src.write_bytes(body)
        out = _apply_file(src, tmp_path / f"{variant}.out.ann",
                          lengths=read_fasta_lengths(FIX / f"{name}.fasta"))
        result = out.read_bytes()
        _assert_feature_rows_equal(src, out)
        assert result.endswith(b"\n") == body.endswith(b"\n"), variant
        if variant.startswith("crlf"):
            assert _no_bare_lf(result), variant


@pytest.mark.parametrize("name", FIXTURES)
def test_roundtrip_meta_equal(name, tmp_path):
    ann, doc, lengths = _load(name)
    out = tmp_path / "out.ann"
    apply(ann, _common(doc), doc.entries, category=_category(doc)).write(out)
    doc2 = extract(AnnFile.read(out), lengths)
    assert doc2.common == doc.common
    assert doc2.common_extra == doc.common_extra
    assert doc2.source == doc.source
    assert doc2.untouched_common_features == doc.untouched_common_features
    if _category(doc) == "GENE":
        assert doc2.entries == doc.entries          # ff_definition is written as-is
    else:
        strip = [dataclasses.replace(e, ff_definition=None) for e in doc.entries]
        assert [dataclasses.replace(e, ff_definition=None) for e in doc2.entries] == strip
    # a second pass changes nothing
    out2 = tmp_path / "out2.ann"
    apply(AnnFile.read(out), _common(doc2), doc2.entries, category=_category(doc)).write(out2)
    assert out2.read_bytes() == out.read_bytes()


# -- layouts and regeneration --------------------------------------------------

def test_dfast_draft_moves_source_to_common(tmp_path):
    ann, doc, _ = _load("dfast_draft")
    assert _category(doc) == "WGS"
    apply(ann, _common(doc), doc.entries, category="WGS")
    feats = ann.features()
    common_src = [f for f in feats if f.entry == "COMMON" and f.key == "source"]
    assert len(common_src) == 1
    src = common_src[0]
    assert src.location == "1..E"
    assert [f for f in feats if f.entry == "COMMON"][-1] is src      # at the end of COMMON
    assert [f for f in feats if f.entry != "COMMON" and f.key == "source"] == []
    assert src.get(ann.rows, "submitter_seqid") == "@@[entry]@@"
    assert src.get(ann.rows, "mol_type") == "genomic DNA"
    # WGS meta definition: the per-contig name comes from submitter_seqid (= @@[entry]@@)
    assert src.get(ann.rows, "ff_definition") == "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"
    # the CDS rows now carry their entry's marker
    assert ann.rows[[f.first for f in feats if f.key == "CDS"][0]].cols[0] == "sequence1"


def test_dfast_draft_with_override_stays_per_entry():
    ann, doc, _ = _load("dfast_draft")
    doc.entries[1].overrides["strain"] = "X-1"
    apply(ann, _common(doc), doc.entries, category="WGS")
    feats = ann.features()
    assert not any(f.entry == "COMMON" and f.key == "source" for f in feats)
    srcs = [f for f in feats if f.key == "source"]
    assert [f.entry for f in srcs] == ["sequence1", "sequence2", "sequence3"]
    assert srcs[0].qualifiers(ann.rows) == [
        ("mol_type", "genomic DNA"), ("submitter_seqid", "sequence1"),
        ("ff_definition", "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@")]
    assert ("strain", "X-1") in srcs[1].qualifiers(ann.rows)


def test_dfast_complete_regenerates_roles():
    ann, doc, _ = _load("dfast_complete")
    doc.entries[1].seq_name = "pNEW"
    apply(ann, _common(doc), doc.entries, category="GNM")
    feats = ann.features()
    by_entry = {e: [f for f in feats if f.entry == e] for e in ("sequence1", "sequence2")}
    for eid, fs in by_entry.items():
        keys = [f.key for f in fs]
        assert keys.count("TOPOLOGY") == 1
        t = keys.index("TOPOLOGY")
        assert keys[t + 1] == "source"                       # TOPOLOGY directly before source
        assert ann.rows[fs[t].first].cols == [eid, "TOPOLOGY", "", "circular", ""]
        assert ann.rows[fs[t + 1].first].cols[0] == ""       # marker moved to TOPOLOGY
    p_src = next(f for f in by_entry["sequence2"] if f.key == "source")
    assert p_src.location == "1..1500"
    assert p_src.get(ann.rows, "plasmid") == "pNEW"
    ff = p_src.get(ann.rows, "ff_definition")
    assert ff == "@@[organism]@@ @@[strain]@@ plasmid @@[plasmid]@@ DNA, complete sequence"
    c_src = next(f for f in by_entry["sequence1"] if f.key == "source")
    quals = c_src.qualifiers(ann.rows)
    names = [q for q, _ in quals]
    # base source (JSON order) -> role -> ff_definition
    assert names == ["organism", "strain", "mol_type", "geo_loc_name", "collection_date", "host",
                     "isolation_source", "chromosome", "ff_definition"]
    assert c_src.get(ann.rows, "ff_definition") == "@@[organism]@@ @@[strain]@@ DNA, chromosome, complete genome"


def test_marker_set_keeps_ff_definition_and_extra():
    ann, doc, _ = _load("marker_set")
    doc.entries[0].ff_definition = "@@[organism]@@ COI gene for cytochrome oxidase, partial cds"
    doc.entries[2].ff_definition = ""
    apply(ann, _common(doc), doc.entries, category="GENE")
    srcs = [f for f in ann.features() if f.key == "source"]
    assert [f.entry for f in srcs] == ["COI_01", "COI_02", "COI_03"]
    q1, q2, q3 = (f.qualifiers(ann.rows) for f in srcs)
    assert ("ff_definition", "@@[organism]@@ COI gene for cytochrome oxidase, partial cds") in q1
    assert not any(q == "ff_definition" for q, _ in q3)
    assert [q for q, _ in q2] == ["mol_type", "organism", "isolate", "collection_date", "geo_loc_name",
                                  "ff_definition", "note", "note"]
    assert [v for q, v in q2 if q == "note"] == ["larva", "voucher deposited in NSMT"]
    assert not any(f.key == "TOPOLOGY" for f in ann.features())


def test_handmade_keeps_unknown_common_feature_and_extra_qualifier(tmp_path):
    ann, doc, _ = _load("handmade")
    out = tmp_path / "out.ann"
    apply(ann, _common(doc), doc.entries, category="WGS").write(out)
    data = out.read_bytes()
    assert b"\n\tKEYWORD\t\tkeyword\tWGS\r\n" in data
    assert _no_bare_lf(data)
    feats = ann.features()
    sub = next(f for f in feats if f.key == "SUBMITTER")
    quals = sub.qualifiers(ann.rows)
    assert quals[-1] == ("fax", "81-55-981-0000")
    assert [q for q, _ in quals[:-1]] == ["ab_name", "ab_name", "contact", "email", "institute",
                                          "country", "city", "street", "zip"]
    assert [f.key for f in feats if f.entry == "COMMON"] == ["SUBMITTER", "REFERENCE", "DATE", "KEYWORD",
                                                             "source"]
    assert ann.rows[0].cols[0] == "COMMON"
    assert sum(1 for r in ann.rows if r.cols[0] == "COMMON") == 1
    src = feats[[f.key for f in feats].index("source")]
    assert src.get(ann.rows, "sub_strain") == "HS-7a"


def test_apply_inserts_dblink_when_absent():
    ann, doc, _ = _load("handmade")
    common = _common(doc)
    common["DBLINK"] = {"project": "PRJDB1", "sample": "SAMD1"}
    apply(ann, common, doc.entries, category="WGS")
    feats = ann.features()
    assert [f.key for f in feats[:3]] == ["DBLINK", "SUBMITTER", "REFERENCE"]
    assert ann.rows[0].cols == ["COMMON", "DBLINK", "", "project", "PRJDB1"]
    assert ann.rows[1].cols == ["", "", "", "biosample", "SAMD1"]
    assert ann.rows[0].eol == "\r\n"
    assert ann.rows[feats[1].first].cols[0] == ""


def test_apply_places_common_targets_after_untouched_leading_features():
    ann, doc, _ = _load("dfast_draft")
    common = _common(doc)
    common["DATATYPE"] = {"type": "WGS"}        # not a target: ignored, the file's own row stays
    apply(ann, common, doc.entries, category="WGS")
    keys = [f.key for f in ann.features() if f.entry == "COMMON"]
    assert keys == ["DATATYPE", "KEYWORD", "DBLINK", "SUBMITTER", "REFERENCE", "COMMENT", "ST_COMMENT",
                    "source"]
    assert sum(1 for r in ann.rows if r.cols[0] == "COMMON") == 1


def _ann(tmp_path, lines: list[str], eol: str = "\n", final: bool = True, name: str = "x.ann") -> Path:
    p = tmp_path / name
    p.write_bytes((eol.join(lines) + (eol if final else "")).encode())
    return p


MINI = [
    "COMMON\tSUBMITTER\t\tab_name\tSuzuki,H.",
    "seq1\tsource\t1..300\torganism\tBacillus subtilis",
    "\t\t\tmol_type\tgenomic DNA",
    "\tCDS\t1..300\tproduct\thypothetical protein",
    "seq2\tCDS\t1..200\tproduct\thypothetical protein",
    "\tsource\t1..200\torganism\tBacillus subtilis",
    "\t\t\tmol_type\tgenomic DNA",
]


def test_apply_keeps_missing_final_newline(tmp_path):
    src = _ann(tmp_path, MINI, final=False)

    def mutate(common, entries):
        entries[1].overrides["strain"] = "HS-9"

    out = _apply_file(src, tmp_path / "out.ann", category="GENE", mutate=mutate)
    data = out.read_bytes()
    assert not data.endswith(b"\n")
    assert data.endswith(b"\t\t\tstrain\tHS-9")
    _assert_feature_rows_equal(src, out)


def test_apply_crlf_without_final_newline_stays_crlf(tmp_path):
    lines = MINI[:5] + ["\tsource\t1..200\torganism\tBacillus subtilis"]   # one-row source at the end
    src = _ann(tmp_path, lines, eol="\r\n", final=False)

    def mutate(common, entries):
        entries[1].overrides["strain"] = "HS-9"
        entries[1].type, entries[1].seq_name, entries[1].status = "plasmid", "pB", "complete"
        entries[1].topology = "circular"

    out = _apply_file(src, tmp_path / "out.ann", category="GNM", mutate=mutate)
    data = out.read_bytes()
    assert not data.endswith(b"\n")
    assert _no_bare_lf(data)
    assert b"\r\n\tTOPOLOGY\t\tcircular\t\r\n\tsource\t1..200\t" in data
    _assert_feature_rows_equal(src, out)


def test_apply_inserts_missing_entry_source(tmp_path):
    lines = MINI[:4] + ["seq2\tCDS\t1..200\tproduct\thypothetical protein"]
    src = _ann(tmp_path, lines)
    ann = AnnFile.read(src)
    doc = extract(ann, {"seq1": 300, "seq2": 200, "seq3": 50})
    doc.entries[2].topology = "circular"
    doc.entries[2].type, doc.entries[2].seq_name, doc.entries[2].status = "plasmid", "p3", "complete"
    apply(ann, _common(doc), doc.entries, category="GNM")
    out = tmp_path / "out.ann"
    ann.write(out)
    feats = ann.features()
    seq2 = [f for f in feats if f.entry == "seq2"]
    assert [f.key for f in seq2] == ["source", "CDS"]
    assert seq2[0].location == "1..200"
    assert ann.rows[seq2[0].first].cols[0] == "seq2" and ann.rows[seq2[1].first].cols[0] == ""
    seq3 = [f for f in feats if f.entry == "seq3"]
    assert [f.key for f in seq3] == ["TOPOLOGY", "source"]
    assert seq3[1].location == "1..50"
    assert seq3[1].get(ann.rows, "plasmid") == "p3"
    _assert_feature_rows_equal(src, out)


def test_common_source_layout_drops_rows_of_featureless_contig(tmp_path):
    lines = (FIX / "dfast_draft.ann").read_text(encoding="utf-8").split("\n")
    s2 = lines.index(next(ln for ln in lines if ln.startswith("sequence2\tsource")))
    s3 = lines.index(next(ln for ln in lines if ln.startswith("sequence3\tsource")))
    lines = lines[:s2 + 5] + lines[s3:]          # sequence2 keeps only its 5 source rows
    src = tmp_path / "in.ann"
    src.write_text("\n".join(lines), encoding="utf-8")
    lengths = read_fasta_lengths(FIX / "dfast_draft.fasta")
    out = _apply_file(src, tmp_path / "out.ann", category="WGS", lengths=lengths)
    feats = AnnFile.read(out).features()
    assert not any(f.entry == "sequence2" for f in feats)
    assert [f.entry for f in feats if f.key == "CDS"] == ["sequence1", "sequence3"]
    _assert_feature_rows_equal(src, out)
    doc = extract(AnnFile.read(out), lengths)
    assert [e.seq_id for e in doc.entries] == ["sequence1", "sequence3", "sequence2"]
    assert doc.layout == "common_source"
    # a second pass keeps the layout and does not append a source for sequence2
    out2 = tmp_path / "out2.ann"
    apply(AnnFile.read(out), _common(doc), doc.entries, category="WGS").write(out2)
    assert out2.read_bytes() == out.read_bytes()


def test_apply_creates_common_block_when_absent(tmp_path):
    ann = AnnFile.read(_ann(tmp_path, MINI[1:]))
    doc = extract(ann)
    common = _common(doc)
    common["SUBMITTER"] = {"ab_name": ["Suzuki,H."], "contact": "Hanako Suzuki"}
    apply(ann, common, doc.entries, category="GENE")
    assert ann.rows[0].cols == ["COMMON", "SUBMITTER", "", "ab_name", "Suzuki,H."]
    assert ann.rows[2].cols[:2] == ["seq1", "source"]


def test_apply_writes_flags_and_drops_false(tmp_path):
    ann = AnnFile.read(_ann(tmp_path, MINI))
    doc = extract(ann)
    common = _common(doc)
    common["SOURCE"]["environmental_sample"] = True
    common["SOURCE"]["focus"] = False
    doc.entries[0].extra.append(("my_flag", True))
    doc.entries[1].overrides["germline"] = "no"
    apply(ann, common, doc.entries, category="GENE")
    s1, s2 = [f for f in ann.features() if f.key == "source"]
    assert ("environmental_sample", "") in s1.qualifiers(ann.rows)
    assert ("my_flag", "") in s1.qualifiers(ann.rows)
    assert not any(q in ("focus", "germline") for q, _ in s1.qualifiers(ann.rows) + s2.qualifiers(ann.rows))


def test_apply_rejects_invalid():
    # The CLI (Task 5) runs validate first and never calls apply on issues.
    _, doc, lengths = _load("handmade")
    common = _common(doc)
    common["SUBMITTER"] = {**common["SUBMITTER"], "contact": "Hanako\tSuzuki"}
    issues = validate(common, doc.entries, category="WGS", fasta_ids=list(lengths),
                      allowed_unknown=set(doc.unknown_qualifiers))
    assert issues == [Issue("SUBMITTER.contact", None, "forbidden_char")]
    common["SUBMITTER"]["contact"] = "Hanako Suzuki"
    assert validate(common, doc.entries, category="WGS", fasta_ids=list(lengths),
                    allowed_unknown=set(doc.unknown_qualifiers)) == []


@pytest.mark.parametrize("name", FIXTURES)
def test_fixtures_validate_clean(name):
    _, doc, lengths = _load(name)
    issues = validate(_common(doc), doc.entries, category=_category(doc), fasta_ids=list(lengths),
                      allowed_unknown=set(doc.unknown_qualifiers))
    assert issues == []


# -- entry TSV -----------------------------------------------------------------

def test_read_entries_tsv_roundtrip(tmp_path):
    entries: list[Entry] = []
    for name in FIXTURES:
        for e in _load(name)[1].entries:
            names = [q for q, _ in e.extra]
            if len(names) != len(set(names)):
                continue                       # repeated extra names cannot be shown in a TSV cell
            entries.append(dataclasses.replace(e, length=None))
    entries.append(Entry("x1", None, type="organelle", seq_name="mitochondrion", status="complete",
                         topology="circular", overrides={"environmental_sample": True, "strain": "A"},
                         ff_definition="@@[organism]@@ mitochondrial DNA, complete genome",
                         extra=[("my_flag", True), ("plasmid", "p1"), ("sub_strain", "s")]))
    out = tmp_path / "e.tsv"
    write_entries(entries, out)
    assert read_entries(out) == entries


def test_read_entries_rejects_bad_tables(tmp_path):
    p = tmp_path / "e.tsv"
    p.write_text("seq1\tunplaced\n", encoding="utf-8")
    with pytest.raises(ValueError, match="header"):
        read_entries(p)
    p.write_text("#seq_id\ttype\tseq_name\tstatus\ttopology\n", encoding="utf-8")
    with pytest.raises(ValueError, match="ff_definition"):
        read_entries(p)
    p.write_text("#seq_id\ttype\tseq_name\tstatus\ttopology\tff_definition\tnote\tnote\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate"):
        read_entries(p)
    p.write_text("#seq_id\ttype\tseq_name\tstatus\ttopology\tff_definition\ns1\tunplaced\t\tpartial\t\t\tx\n",
                 encoding="utf-8")
    with pytest.raises(ValueError, match="line 2"):
        read_entries(p)


def test_read_entries_crlf_short_rows_and_blank_lines(tmp_path):
    p = tmp_path / "e.tsv"
    p.write_bytes(b"#seq_id\ttype\tseq_name\tstatus\ttopology\tff_definition\tstrain\r\n"
                  b"s1\tchromosome\t1\tcomplete\tcircular\r\n\r\ns2\r\n")
    assert read_entries(p) == [
        Entry("s1", None, type="chromosome", seq_name="1", status="complete", topology="circular"),
        Entry("s2", None)]


# -- fix round 1 -----------------------------------------------------------------

def test_apply_raises_for_ann_entry_missing_from_table():
    for name, category in (("dfast_draft", "WGS"), ("dfast_complete", "GNM")):
        ann, doc, _ = _load(name)
        before = [list(r.cols) + [r.eol] for r in ann.rows]
        with pytest.raises(ValueError, match="entry_missing_from_table: sequence2"):
            apply(ann, _common(doc), [e for e in doc.entries if e.seq_id != "sequence2"], category=category)
        assert [list(r.cols) + [r.eol] for r in ann.rows] == before      # nothing touched


def test_common_source_layout_removes_every_entry_source(tmp_path):
    lines = [
        "COMMON\tsource\t1..E\torganism\tBacillus subtilis",
        "\t\t\tmol_type\tgenomic DNA",
        "seq1\tsource\t1..300\torganism\tBacillus subtilis",
        "\t\t\tmol_type\tgenomic DNA",
        "\tTOPOLOGY\t\tlinear\t",
        "\tCDS\t1..300\tproduct\thypothetical protein",
        "seq2\tCDS\t1..200\tproduct\thypothetical protein",
    ]
    src = _ann(tmp_path, lines)
    lengths = {"seq1": 300, "seq2": 200}
    out = _apply_file(src, tmp_path / "out.ann", category="WGS", lengths=lengths)
    feats = AnnFile.read(out).features()
    assert [(f.entry, f.key) for f in feats] == [("COMMON", "source"), ("seq1", "CDS"), ("seq2", "CDS")]
    doc = extract(AnnFile.read(out), lengths)
    assert not any(n.code == "both_common_and_entry_source" for n in doc.notes)
    out2 = tmp_path / "out2.ann"
    apply(AnnFile.read(out), _common(doc), doc.entries, category="WGS").write(out2)
    assert out2.read_bytes() == out.read_bytes()


def test_literal_yes_survives_tsv_and_apply(tmp_path):
    ann = AnnFile.read(_ann(tmp_path, MINI))
    doc = extract(ann)
    doc.entries[0].overrides["strain"] = "yes"
    doc.entries[0].overrides["environmental_sample"] = True
    doc.entries[0].extra.append(("my_flag", True))
    tsv = tmp_path / "e.tsv"
    write_entries(doc.entries, tsv)
    entries = read_entries(tsv)
    assert entries[0].overrides["strain"] == "yes"
    assert entries[0].overrides["environmental_sample"] is True
    assert entries[0].extra == [("my_flag", True)]
    for e, orig in zip(entries, doc.entries):
        e.length = orig.length
    apply(ann, _common(doc), entries, category="GENE")
    s1 = next(f for f in ann.features() if f.key == "source")
    quals = s1.qualifiers(ann.rows)
    assert ("strain", "yes") in quals
    assert ("environmental_sample", "") in quals and ("my_flag", "") in quals


# -- final review fixes ------------------------------------------------------------

TSA_ANN = [
    "COMMON\tDATATYPE\t\ttype\tTSA",
    "\tsource\t1..E\torganism\tOryza sativa",
    "\t\t\tmol_type\ttranscribed RNA",
    "\t\t\tsubmitter_seqid\t@@[entry]@@",
    "\t\t\tff_definition\t@@[organism]@@ RNA, @@[submitter_seqid]@@",
    "t1\tCDS\t1..30\tproduct\tp",
]


@pytest.mark.parametrize("category", ["TSA", "TSA-ENV"])
def test_tsa_uses_the_wgs_layout(tmp_path, category):
    # I4: TSA follows tsa2mss -- one COMMON source, contigs named by submitter_seqid
    ann = AnnFile.read(_ann(tmp_path, TSA_ANN))
    doc = extract(ann, {"t1": 300, "t2": 300})
    apply(ann, _common(doc), doc.entries, category=category)
    feats = ann.features()
    assert [(f.entry, f.key) for f in feats] == [("COMMON", "DATATYPE"), ("COMMON", "source"), ("t1", "CDS")]
    quals = feats[1].qualifiers(ann.rows)
    assert ("submitter_seqid", "@@[entry]@@") in quals
    assert ("ff_definition", "@@[organism]@@ RNA, @@[submitter_seqid]@@") in quals
    assert (("environmental_sample", "") in quals) == (category == "TSA-ENV")


def test_tsa_per_entry_names_contigs_by_submitter_seqid(tmp_path):
    ann = AnnFile.read(_ann(tmp_path, TSA_ANN))
    doc = extract(ann, {"t1": 300, "t2": 300})
    doc.entries[1].overrides["note"] = "x"
    apply(ann, _common(doc), doc.entries, category="TSA")
    srcs = [f for f in ann.features() if f.key == "source"]
    assert [f.entry for f in srcs] == ["t1", "t2"]
    q1 = srcs[0].qualifiers(ann.rows)
    assert ("submitter_seqid", "t1") in q1
    assert ("ff_definition", "@@[organism]@@ RNA, @@[submitter_seqid]@@") in q1
