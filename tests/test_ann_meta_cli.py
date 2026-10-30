import json
import shutil
from pathlib import Path

import pytest

from ann_fix.annfile import AnnFile
from ann_meta.cli import main
from ann_meta.entries_tsv import read_entries
from ann_meta.extract import COMMON_TARGETS

FIX = Path(__file__).parent / "ann_meta_fixtures"


def _common_from_doc(doc: dict) -> dict:
    """The common JSON as GGS assembles it from `ann_meta extract` output."""
    c = dict(doc["common"])
    c["SOURCE"] = dict(doc["source"])
    if doc["source_identifier"]:
        c["SOURCE_IDENTIFIER"] = doc["source_identifier"]
    if doc["common_extra"]:
        c["_extra"] = doc["common_extra"]
    return c


def _extract(tmp_path, capsys, name="dfast_complete"):
    tsv = tmp_path / "entries.tsv"
    rc = main(["extract", "--ann", str(FIX / f"{name}.ann"), "--fasta", str(FIX / f"{name}.fasta"),
               "--entries-out", str(tsv)])
    out = capsys.readouterr().out
    assert rc == 0
    doc = json.loads(out)
    common = tmp_path / "common.json"
    common.write_text(json.dumps(_common_from_doc(doc)), encoding="utf-8")
    return doc, common, tsv


def _feature_rows(path: Path) -> list[bytes]:
    out, entry, key = [], "", ""
    for raw in path.read_bytes().splitlines():
        cols = raw.split(b"\t") + [b""] * 5
        if cols[0]:
            entry = cols[0].decode()
        if cols[1]:
            key = cols[1].decode()
        if key in ("source", "TOPOLOGY") or (entry == "COMMON" and key in COMMON_TARGETS):
            continue
        out.append(raw)
    return out


def test_extract_prints_json_and_writes_tsv(tmp_path, capsys):
    doc, _, tsv = _extract(tmp_path, capsys)
    assert doc["category"]["value"] == "GNM"
    assert doc["common"]["DBLINK"]["project"] == "PRJDB99901"
    ids = [e["seq_id"] for e in doc["entries"]]
    assert ids and ids[0] == "sequence1"
    assert [e.seq_id for e in read_entries(tsv)] == ids
    # --json is accepted for symmetry; output is the same JSON
    assert main(["extract", "--ann", str(FIX / "dfast_complete.ann"), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["entries"][0]["seq_id"] == "sequence1"


def test_apply_roundtrip_exit_0(tmp_path, capsys):
    _, common, tsv = _extract(tmp_path, capsys)
    out = tmp_path / "out.ann"
    rc = main(["apply", "--ann", str(FIX / "dfast_complete.ann"), "--common", str(common),
               "--entries", str(tsv), "--out", str(out), "--category", "GNM",
               "--fasta", str(FIX / "dfast_complete.fasta")])
    assert rc == 0, capsys.readouterr()
    assert _feature_rows(out) == _feature_rows(FIX / "dfast_complete.ann")
    # without --fasta, lengths come from the existing sources; still a clean round trip
    out2 = tmp_path / "out2.ann"
    assert main(["apply", "--ann", str(FIX / "dfast_complete.ann"), "--common", str(common),
                 "--entries", str(tsv), "--out", str(out2), "--category", "GNM"]) == 0
    assert out2.read_bytes() == out.read_bytes()


def test_apply_validation_issue_exit_2_and_leaves_output_untouched(tmp_path, capsys):
    _, common, tsv = _extract(tmp_path, capsys)
    data = json.loads(common.read_text())
    data["SUBMITTER"]["contact"] = "Taro\tMishima"
    common.write_text(json.dumps(data), encoding="utf-8")
    out = tmp_path / "out.ann"
    out.write_bytes(b"untouched")
    rc = main(["apply", "--ann", str(FIX / "dfast_complete.ann"), "--common", str(common),
               "--entries", str(tsv), "--out", str(out), "--category", "GNM"])
    assert rc == 2
    res = json.loads(capsys.readouterr().out)
    assert {"field": "SUBMITTER.contact", "entry": None, "message": "forbidden_char"} in res["issues"]
    assert out.read_bytes() == b"untouched"


def test_apply_in_place(tmp_path, capsys):
    _, common, tsv = _extract(tmp_path, capsys)
    ann = tmp_path / "x.ann"
    shutil.copy(FIX / "dfast_complete.ann", ann)
    data = json.loads(common.read_text())
    data["SUBMITTER"]["contact"] = "Hanako Mishima"
    common.write_text(json.dumps(data), encoding="utf-8")
    rc = main(["apply", "--ann", str(ann), "--common", str(common), "--entries", str(tsv),
               "--out", str(ann), "--category", "GNM"])
    assert rc == 0
    rows = AnnFile.read(ann).rows
    assert any(r.cols[3:5] == ["contact", "Hanako Mishima"] for r in rows)
    assert _feature_rows(ann) == _feature_rows(FIX / "dfast_complete.ann")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["common.json", "entries.tsv", "x.ann"]


def test_extract_unreadable_exit_1(tmp_path, capsys):
    rc = main(["extract", "--ann", str(tmp_path / "missing.ann")])
    assert rc == 1
    err = capsys.readouterr().err
    assert err.startswith("error:") and "Traceback" not in err
    bad = tmp_path / "bin.ann"
    bad.write_bytes(b"\xff\xfe\x00\x01")
    assert main(["extract", "--ann", str(bad)]) == 1


# -- extra rejection cases (controller rulings) ----------------------------------

def _apply(tmp_path, common, tsv, *extra):
    return main(["apply", "--ann", str(FIX / "dfast_complete.ann"), "--common", str(common),
                 "--entries", str(tsv), "--out", str(tmp_path / "o.ann"), "--category", "GNM", *extra])


def test_apply_bad_json_exit_1(tmp_path, capsys):
    _, _, tsv = _extract(tmp_path, capsys)
    bad = tmp_path / "bad.json"
    bad.write_text("[1, 2]", encoding="utf-8")
    assert _apply(tmp_path, bad, tsv) == 1
    bad.write_text("{not json", encoding="utf-8")
    assert _apply(tmp_path, bad, tsv) == 1
    assert _apply(tmp_path, tmp_path / "nope.json", tsv) == 1
    assert not (tmp_path / "o.ann").exists()


@pytest.mark.parametrize("mutate", [
    lambda d: d.__setitem__("REFERENCE", "not a list"),
    lambda d: d.__setitem__("_extra", {"COMMENT": [["line", 3]]}),
    lambda d: d.__setitem__("_extra", ["x"]),
    lambda d: d.__setitem__("_extra", {"COMMENT": [["a", "b", "c"]]}),
])
def test_apply_bad_shape_exit_2(tmp_path, capsys, mutate):
    _, common, tsv = _extract(tmp_path, capsys)
    data = json.loads(common.read_text())
    mutate(data)
    common.write_text(json.dumps(data), encoding="utf-8")
    assert _apply(tmp_path, common, tsv) == 2
    issues = json.loads(capsys.readouterr().out)["issues"]
    assert issues and all(i["message"].startswith("bad_shape:") and i["entry"] is None for i in issues)


def test_apply_bad_entries_tsv_exit_2(tmp_path, capsys):
    _, common, _ = _extract(tmp_path, capsys)
    tsv = tmp_path / "bad.tsv"
    tsv.write_text("seq_id\ttype\n", encoding="utf-8")
    assert _apply(tmp_path, common, tsv) == 2
    issues = json.loads(capsys.readouterr().out)["issues"]
    assert issues[0]["message"].startswith("bad_entries_tsv")


def test_apply_entry_missing_from_table_exit_2(tmp_path, capsys):
    _, common, tsv = _extract(tmp_path, capsys)
    lines = tsv.read_text().splitlines()
    tsv.write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")   # drop the last entry
    assert _apply(tmp_path, common, tsv) == 2
    issues = json.loads(capsys.readouterr().out)["issues"]
    assert any(i["message"].startswith("entry_missing_from_table") for i in issues)
    assert not (tmp_path / "o.ann").exists()


def test_apply_allow_unknown_explicit(tmp_path, capsys):
    _, common, tsv = _extract(tmp_path, capsys)
    lines = tsv.read_text().splitlines()
    lines[0] += "\tmy_custom_q"
    lines[1] += "\tv1"
    tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")
    assert _apply(tmp_path, common, tsv) == 2
    assert {"field": "my_custom_q", "entry": "sequence1", "message": "unknown_qualifier"} in \
        json.loads(capsys.readouterr().out)["issues"]
    assert _apply(tmp_path, common, tsv, "--allow-unknown", "my_custom_q") == 0


def test_apply_extra_flag_values(tmp_path, capsys):
    _, common, tsv = _extract(tmp_path, capsys)
    data = json.loads(common.read_text())
    data["_extra"] = {"COMMENT": [["flag_on", True], ["flag_off", False]]}
    common.write_text(json.dumps(data), encoding="utf-8")
    assert _apply(tmp_path, common, tsv) == 0
    rows = [r.cols[3:5] for r in AnnFile.read(tmp_path / "o.ann").rows]
    assert ["flag_on", ""] in rows
    assert not any(c[0] == "flag_off" for c in rows)


def test_apply_rejects_extra_for_non_target_feature(tmp_path, capsys):
    _, common, tsv = _extract(tmp_path, capsys)
    data = json.loads(common.read_text())
    data["_extra"] = {"KEYWORD": [["keyword", "x"]], "COMMENT": [["line", "ok"]]}
    common.write_text(json.dumps(data), encoding="utf-8")
    assert _apply(tmp_path, common, tsv) == 2
    issues = json.loads(capsys.readouterr().out)["issues"]
    assert issues == [{"field": "_extra.KEYWORD", "entry": None, "message": "bad_shape: not a COMMON target"}]
    assert not (tmp_path / "o.ann").exists()


# -- fix round 1 -------------------------------------------------------------------

def test_apply_in_place_keeps_file_mode(tmp_path, capsys):
    _, common, tsv = _extract(tmp_path, capsys)
    ann = tmp_path / "x.ann"
    shutil.copy(FIX / "dfast_complete.ann", ann)
    ann.chmod(0o644)
    assert main(["apply", "--ann", str(ann), "--common", str(common), "--entries", str(tsv),
                 "--out", str(ann), "--category", "GNM"]) == 0
    assert ann.stat().st_mode & 0o7777 == 0o644


def test_apply_new_out_gets_umask_default(tmp_path, capsys):
    import os
    _, common, tsv = _extract(tmp_path, capsys)
    old = os.umask(0o027)
    try:
        assert _apply(tmp_path, common, tsv) == 0
    finally:
        os.umask(old)
    assert (tmp_path / "o.ann").stat().st_mode & 0o7777 == 0o640


def test_usage_error_exit_1_empty_stdout(tmp_path, capsys):
    rc = main(["apply", "--ann", "x.ann", "--common", "c.json", "--entries", "e.tsv", "--out", "o.ann"])
    assert rc == 1
    cap = capsys.readouterr()
    assert cap.out == ""
    assert "--category" in cap.err


def test_internal_error_exit_1_no_traceback(tmp_path, capsys, monkeypatch):
    import ann_meta.cli as cli
    _, common, tsv = _extract(tmp_path, capsys)

    def boom(*a, **k):
        raise IndexError("list index out of range")
    monkeypatch.setattr(cli, "apply", boom)
    assert _apply(tmp_path, common, tsv) == 1
    cap = capsys.readouterr()
    assert cap.out == ""
    assert cap.err.strip() == "error: internal: IndexError: list index out of range"
    assert not (tmp_path / "o.ann").exists()


# -- final review fixes ------------------------------------------------------------

def _mutated(tmp_path, capsys, mutate):
    _, common, tsv = _extract(tmp_path, capsys)
    data = json.loads(common.read_text())
    mutate(data)
    common.write_text(json.dumps(data), encoding="utf-8")
    return common, tsv


def _rows(path: Path) -> list[list[str]]:
    return [r.cols[3:5] for r in AnnFile.read(path).rows]


@pytest.mark.parametrize("mutate", [
    lambda d: d.__setitem__("DBLINK", {"project": "PRJDB1"}),
    lambda d: d["SUBMITTER"].pop("ab_name"),
    lambda d: d["SUBMITTER"].__setitem__("ab_name", []),
    lambda d: d["REFERENCE"][0].update(year="", status="unpublished"),
    lambda d: d["REFERENCE"][0].update(status="Published", journal="J"),
    lambda d: d["REFERENCE"][0].pop("title"),
    lambda d: d.__setitem__("SUBMITTER", {"contact": "X", "fax": "1"}),
], ids=["dblink_project_only", "submitter_no_ab_name", "submitter_empty_ab_name",
        "reference_blank_year_lower_status", "published_without_volume", "reference_no_title",
        "submitter_extra_field"])
def test_apply_accepts_partially_filled_common(tmp_path, capsys, mutate):
    # I2: required items are not a gate; the shape check is structural only
    common, tsv = _mutated(tmp_path, capsys, mutate)
    assert _apply(tmp_path, common, tsv) == 0, capsys.readouterr().out


def test_apply_drops_null_fields(tmp_path, capsys):
    common, tsv = _mutated(tmp_path, capsys, lambda d: d["SUBMITTER"].update(contact=None, fax=None))
    assert _apply(tmp_path, common, tsv) == 0, capsys.readouterr().out
    rows = _rows(tmp_path / "o.ann")
    assert not any(q in ("contact", "fax") for q, _ in rows)
    assert not any(v == "None" for _, v in rows)


@pytest.mark.parametrize("mutate", [
    lambda d: d["REFERENCE"][0].__setitem__("year", 2026),
    lambda d: d["SUBMITTER"].__setitem__("contact", 5),
    lambda d: d["SUBMITTER"].__setitem__("ab_name", "Mishima,T."),
    lambda d: d["SUBMITTER"].__setitem__("fax", {"a": "b"}),
    lambda d: d.__setitem__("DBLINK", ["PRJDB1"]),
    lambda d: d["SOURCE"].__setitem__("note", 5),
    lambda d: d.__setitem__("SOURCE_IDENTIFIER", ["strain"]),
    lambda d: d.__setitem__("_extra", {"COMMENT": [["", "x"]]}),
    lambda d: d.__setitem__("DBLINK", {"project": "PRJDB1", "bioproject": "PRJDB2"}),
    lambda d: d.__setitem__("DBLINK", {"sra": ["DRR1"], "sequence read archive": ["DRR2"]}),
], ids=["year_int", "contact_int", "ab_name_str", "extra_field_dict", "dblink_list", "source_int",
        "source_identifier_list", "extra_empty_name", "dblink_alias_clash", "dblink_sra_alias_clash"])
def test_apply_bad_shape_structural(tmp_path, capsys, mutate):
    common, tsv = _mutated(tmp_path, capsys, mutate)
    assert _apply(tmp_path, common, tsv) == 2
    issues = json.loads(capsys.readouterr().out)["issues"]
    assert issues and all(i["message"].startswith("bad_shape:") for i in issues), issues
    assert not (tmp_path / "o.ann").exists()


def test_apply_maps_dblink_aliases(tmp_path, capsys):
    # M1: COMMON aliases are written with the MSS names
    common, tsv = _mutated(tmp_path, capsys, lambda d: d.__setitem__(
        "DBLINK", {"bioproject": "PRJDB9", "biosample": "SAMD9", "sra": ["DRR9"]}))
    assert _apply(tmp_path, common, tsv) == 0, capsys.readouterr().out
    rows = AnnFile.read(tmp_path / "o.ann").rows
    assert [r.cols for r in rows[:3]] == [["COMMON", "DBLINK", "", "project", "PRJDB9"],
                                          ["", "", "", "biosample", "SAMD9"],
                                          ["", "", "", "sequence read archive", "DRR9"]]


def test_apply_lone_surrogate_exit_2(tmp_path, capsys):
    # M2: a JSON "\ud800" escape is an input issue, not an internal error
    _, common, tsv = _extract(tmp_path, capsys)
    data = json.loads(common.read_text())
    data["SUBMITTER"]["contact"] = "\ud800"
    data["SUBMITTER"]["x\udfff"] = "v"
    common.write_text(json.dumps(data), encoding="utf-8")        # ensure_ascii: stays an escape
    assert _apply(tmp_path, common, tsv) == 2
    issues = json.loads(capsys.readouterr().out)["issues"]
    assert {"field": "SUBMITTER.contact", "entry": None, "message": "forbidden_char"} in issues
    assert any(i["field"].startswith("SUBMITTER.x") and i["message"] == "forbidden_char" for i in issues)
    assert not (tmp_path / "o.ann").exists()


def test_apply_bad_source_identifier_exit_2(tmp_path, capsys):
    # I6
    common, tsv = _mutated(tmp_path, capsys, lambda d: d.__setitem__("SOURCE_IDENTIFIER", "]@@ x @@["))
    assert _apply(tmp_path, common, tsv) == 2
    assert json.loads(capsys.readouterr().out)["issues"] == [
        {"field": "SOURCE_IDENTIFIER", "entry": None, "message": "bad_source_identifier"}]


def test_apply_role_column_in_tsv_is_reserved(tmp_path, capsys):
    # I5: /plasmid via a TSV column on a WGS contig
    tsv = tmp_path / "e.tsv"
    assert main(["extract", "--ann", str(FIX / "dfast_draft.ann"), "--entries-out", str(tsv)]) == 0
    doc = json.loads(capsys.readouterr().out)
    common = tmp_path / "c.json"
    common.write_text(json.dumps(_common_from_doc(doc)), encoding="utf-8")
    lines = tsv.read_text().splitlines()
    lines[0] += "\tplasmid"
    lines[1] += "\tpEVIL"
    tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")
    rc = main(["apply", "--ann", str(FIX / "dfast_draft.ann"), "--common", str(common), "--entries", str(tsv),
               "--out", str(tmp_path / "o.ann"), "--category", "WGS"])
    assert rc == 2
    assert json.loads(capsys.readouterr().out)["issues"] == [
        {"field": "plasmid", "entry": "sequence1", "message": "reserved_qualifier"}]
