from pathlib import Path

from ann_fix.annfile import AnnFile

FX = Path(__file__).parent / "ann_fix_fixtures" / "mini.ann"


def _read_write(tmp_path, data: bytes) -> bytes:
    p = tmp_path / "in.ann"
    p.write_bytes(data)
    AnnFile.read(p).write(tmp_path / "out.ann")
    return (tmp_path / "out.ann").read_bytes()


def test_roundtrip_identical(tmp_path):
    data = FX.read_bytes()
    assert _read_write(tmp_path, data) == data
    crlf = data.replace(b"\n", b"\r\n")
    assert _read_write(tmp_path, crlf) == crlf
    nofinal = data.rstrip(b"\n")
    assert _read_write(tmp_path, nofinal) == nofinal


def test_features_and_line_lookup():
    ann = AnnFile.read(FX)
    feats = [f for f in ann.features() if f.entry != "COMMON"]
    assert [f.key for f in feats] == ["source", "mRNA", "CDS", "mRNA", "CDS", "tRNA"]
    assert all(f.entry == "NC_001133.9" for f in feats)
    cds2 = feats[4]
    assert cds2.get(ann.rows, "locus_tag") == "SCEX_000020"
    assert cds2.location == "538..792"
    assert ann.feature_at_line(cds2.first + 1) is not None
    assert ann.feature_at_line(cds2.first + 1).first == cds2.first
    assert ann.feature_at_line(cds2.first + 2) is None          # a qualifier row is not a feature start
    assert [k for k, _ in cds2.qualifiers(ann.rows)][:3] == ["locus_tag", "transl_table", "codon_start"]
    common = [f for f in ann.features() if f.entry == "COMMON"]
    assert common and common[0].key == "DBLINK"


def test_delete_first_feature_moves_entry_id():
    ann = AnnFile.read(FX)
    source = next(f for f in ann.features() if f.key == "source")
    ann.delete_feature(source)
    first = next(f for f in ann.features() if f.entry == "NC_001133.9")
    assert first.key == "mRNA" and ann.rows[first.first].cols[0] == "NC_001133.9"
    assert not any(f.key == "source" for f in ann.features())


def test_replace_and_add_qualifier():
    ann = AnnFile.read(FX)
    n = len(ann.rows)
    cds2 = [f for f in ann.features() if f.key == "CDS"][1]
    ann.replace_feature(cds2, "misc_feature", cds2.location, [("locus_tag", "SCEX_000020"), ("note", "x")])
    assert len(ann.rows) == n - 3                                  # 5 qualifier rows -> 2
    misc = next(f for f in ann.features() if f.key == "misc_feature")
    assert misc.location == "538..792" and misc.qualifiers(ann.rows) == [("locus_tag", "SCEX_000020"), ("note", "x")]
    trna = next(f for f in ann.features() if f.key == "tRNA")
    ann.add_qualifier(trna, "artificial_location", "low-quality sequence region")
    assert ann.rows[trna.last + 1].cols == ["", "", "", "artificial_location", "low-quality sequence region"]
    assert ann.rows[trna.last + 1].eol == "\n"


def _serialize(ann: AnnFile) -> bytes:
    return "".join("\t".join(r.cols) + r.eol for r in ann.rows).encode("utf-8")


def test_insert_feature_middle(tmp_path):
    data = FX.read_bytes()
    ann = AnnFile.read(FX)
    first = ann.features()[0]
    feat = ann.insert_feature(first.first, "COMMON", "DBLINK", "",
                              [("project", "PRJDB1"), ("biosample", "SAMD1")])
    assert ann.features()[0].key == "DBLINK"
    assert (feat.key, feat.first, feat.last) == ("DBLINK", 0, 1)
    assert ann.rows[0].cols == ["COMMON", "DBLINK", "", "project", "PRJDB1"]
    assert ann.rows[1].cols == ["", "", "", "biosample", "SAMD1"]
    assert ann.rows[0].eol == ann.rows[1].eol == "\n"
    ann.write(tmp_path / "out.ann")
    out = (tmp_path / "out.ann").read_bytes()
    lines = out.splitlines(keepends=True)
    orig = data.splitlines(keepends=True)
    assert lines[2] == orig[0].replace(b"COMMON", b"", 1)       # marker moved to the inserted block
    assert b"".join(lines[3:]) == b"".join(orig[1:])            # every other row unchanged


def test_insert_feature_at_end_keeps_missing_newline(tmp_path):
    nofinal = FX.read_bytes().rstrip(b"\n")
    p = tmp_path / "in.ann"
    p.write_bytes(nofinal)
    ann = AnnFile.read(p)
    n = len(ann.rows)
    feat = ann.insert_feature(n, "", "misc_feature", "1..10", [("note", "a"), ("note", "b")])
    assert (feat.first, feat.last) == (n, n + 1)
    ann.write(tmp_path / "out.ann")
    out = (tmp_path / "out.ann").read_bytes()
    assert not out.endswith(b"\n")
    assert out == nofinal + b"\n\tmisc_feature\t1..10\tnote\ta\n\t\t\tnote\tb"


def test_insert_feature_crlf(tmp_path):
    crlf = FX.read_bytes().replace(b"\n", b"\r\n")
    p = tmp_path / "in.ann"
    p.write_bytes(crlf)
    ann = AnnFile.read(p)
    trna = next(f for f in ann.features() if f.key == "tRNA")
    ann.insert_feature(trna.first, "", "gene", "1..10", [])
    ann.insert_feature(len(ann.rows), "", "misc_feature", "5..6", [("note", "x")])
    ann.write(tmp_path / "out.ann")
    out = (tmp_path / "out.ann").read_bytes()
    assert b"\tgene\t1..10\t\t\r\n" in out
    assert out.endswith(b"\tmisc_feature\t5..6\tnote\tx\r\n")
    assert b"\n" not in out.replace(b"\r\n", b"")
    gene = next(f for f in ann.features() if f.key == "gene")
    assert gene.qualifiers(ann.rows) == []
    assert ann.rows[gene.last + 1].cols[1] == "tRNA"


def test_insert_feature_at_start_keeps_crlf(tmp_path):
    crlf = FX.read_bytes().replace(b"\n", b"\r\n")
    p = tmp_path / "in.ann"
    p.write_bytes(crlf)
    ann = AnnFile.read(p)
    ann.insert_feature(0, "COMMON", "DBLINK", "", [("project", "PRJDB1"), ("biosample", "SAMD1")])
    assert ann.rows[0].eol == ann.rows[1].eol == "\r\n"
    ann.write(tmp_path / "out.ann")
    out = (tmp_path / "out.ann").read_bytes()
    assert out.startswith(b"COMMON\tDBLINK\t\tproject\tPRJDB1\r\n\t\t\tbiosample\tSAMD1\r\n")
    assert b"\n" not in out.replace(b"\r\n", b"")


def test_insert_feature_before_entry_first_row_moves_marker():
    ann = AnnFile.read(FX)
    common_first = next(i for i, r in enumerate(ann.rows) if r.cols[0] == "COMMON")
    feat = ann.insert_feature(common_first, "COMMON", "DBLINK", "", [("project", "PRJDB1")])
    marked = [i for i, r in enumerate(ann.rows) if r.cols[0] == "COMMON"]
    assert marked == [feat.first]
    feats = ann.features()
    assert (feats[0].key, feats[0].entry, feats[0].first) == ("DBLINK", "COMMON", feat.first)
    assert feats[1].entry == "COMMON" and ann.rows[feats[1].first].cols[0] == ""


def test_insert_feature_before_other_entry_keeps_its_marker():
    ann = AnnFile.read(FX)
    src = next(f for f in ann.features() if f.key == "source")
    assert ann.rows[src.first].cols[0] == "NC_001133.9"
    feat = ann.insert_feature(src.first, "", "misc_feature", "1..2", [("note", "x")])
    assert ann.rows[feat.last + 1].cols[0] == "NC_001133.9"
    assert ann.rows[feat.first].cols[0] == ""
