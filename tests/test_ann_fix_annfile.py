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
