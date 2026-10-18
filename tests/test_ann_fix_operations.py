from pathlib import Path

import pytest

from ann_fix.annfile import AnnFile
from ann_fix.operations import (ARTIFICIAL_LOCATION_DEFAULT, DEFAULTS, IGNORABLE, apply, op_artificial_location,
                                op_drop, op_misc, op_overlap_drop, op_overlap_misc, paired_mrna)
from ann_fix.report import Finding

FX = Path(__file__).parent / "ann_fix_fixtures" / "mini.ann"


def _ann():
    return AnnFile.read(FX)


def _keys(ann):
    return [f.key for f in ann.features() if f.entry != "COMMON"]


def _cds2(ann):
    return [f for f in ann.features() if f.key == "CDS"][1]


def _trna(ann):
    return next(f for f in ann.features() if f.key == "tRNA")


def test_defaults_table():
    assert DEFAULTS["ANN5320"] == ("overlap-drop", {"keep": "cds"})
    assert DEFAULTS["ANN4300"] == ("misc", {})
    assert DEFAULTS["AXS6810"] == ("artificial-location", {"value": "low-quality sequence region"})
    assert set(IGNORABLE) == {"ANN2520", "ANN2542"}


def test_paired_mrna_uses_submitter_transcript_id():
    ann = _ann()
    m = paired_mrna(ann, _cds2(ann))
    assert m is not None and m.key == "mRNA" and m.get(ann.rows, "locus_tag") == "SCEX_000020"
    assert paired_mrna(ann, _trna(ann)) is None


def test_misc_converts_cds_and_drops_mrna():
    ann = _ann()
    op_misc(ann, _cds2(ann), reason="ANN4300: Complete CDS must be at least 90 bases")
    assert _keys(ann) == ["source", "mRNA", "CDS", "misc_feature", "tRNA"]
    misc = next(f for f in ann.features() if f.key == "misc_feature")
    quals = misc.qualifiers(ann.rows)
    assert misc.location == "538..792"
    assert quals[0] == ("locus_tag", "SCEX_000020")
    assert quals[1] == ("note", "nonfunctional CDS (ANN4300: Complete CDS must be at least 90 bases)")
    assert quals[2][0] == "note" and quals[2][1].startswith("submitter_gene_id")
    assert [k for k, _ in quals] == ["locus_tag", "note", "note"]


def test_misc_without_mrna():
    ann = _ann()
    op_drop(ann, paired_mrna(ann, _cds2(ann)), pair=False)
    op_misc(ann, _cds2(ann), reason="x")
    assert _keys(ann) == ["source", "mRNA", "CDS", "misc_feature", "tRNA"]


def test_artificial_location_on_cds_and_mrna_idempotent():
    ann = _ann()
    op_artificial_location(ann, _cds2(ann))
    op_artificial_location(ann, _cds2(ann))
    cds = _cds2(ann)
    mrna = paired_mrna(ann, cds)
    assert cds.qualifiers(ann.rows).count(("artificial_location", ARTIFICIAL_LOCATION_DEFAULT)) == 1
    assert mrna.qualifiers(ann.rows).count(("artificial_location", ARTIFICIAL_LOCATION_DEFAULT)) == 1
    assert _keys(ann) == ["source", "mRNA", "CDS", "mRNA", "CDS", "tRNA"]


def test_overlap_drop_keep_cds_removes_trna():
    ann = _ann()
    f = Finding("ANN5320", "ERR", "NC_001133.9", "tRNA", _trna(ann).first + 1,
                "The tRNA feature must not be completely contained within CDS exons. (Found: tRNA 600..671 in CDS 538..792)")
    op_overlap_drop(ann, _trna(ann), finding=f, keep="cds")
    assert _keys(ann) == ["source", "mRNA", "CDS", "mRNA", "CDS"]


def test_overlap_drop_keep_other_removes_cds_pair():
    ann = _ann()
    f = Finding("ANN5320", "ERR", "NC_001133.9", "tRNA", _trna(ann).first + 1,
                "The tRNA feature must not be completely contained within CDS exons. (Found: tRNA 600..671 in CDS 538..792)")
    op_overlap_drop(ann, _trna(ann), finding=f, keep="other")
    assert _keys(ann) == ["source", "mRNA", "CDS", "tRNA"]


def test_overlap_misc_keep_cds_demotes_trna():
    ann = _ann()
    f = Finding("ANN5320", "ERR", "NC_001133.9", "tRNA", _trna(ann).first + 1, "... (Found: tRNA 600..671 in CDS 538..792)")
    op_overlap_misc(ann, _trna(ann), finding=f, keep="cds")
    assert _keys(ann) == ["source", "mRNA", "CDS", "mRNA", "CDS", "misc_feature"]
    misc = next(f for f in ann.features() if f.key == "misc_feature")
    q = dict(misc.qualifiers(ann.rows))
    assert q["locus_tag"] == "SCEX_000030" and "product" not in q and "tRNA-Pro" in " ".join(v for k, v in misc.qualifiers(ann.rows) if k == "note")


def test_apply_descending_and_skips_stale():
    ann = _ann()
    t = _trna(ann).first + 1
    c = _cds2(ann).first + 1
    findings = [Finding("ANN4300", "WAR", "NC_001133.9", "CDS", c, "Complete CDS must be at least 90 bases"),
                Finding("ANN5320", "ERR", "NC_001133.9", "tRNA", t, "(Found: tRNA 600..671 in CDS 538..792)"),
                Finding("XXX9999", "WAR", "NC_001133.9", "tRNA", t, "fake")]
    log = apply(ann, findings, {"ANN4300": ("misc", {}), "ANN5320": ("overlap-drop", {"keep": "cds"}),
                                "XXX9999": ("drop", {})})
    assert _keys(ann) == ["source", "mRNA", "CDS", "misc_feature"]
    assert any(l.startswith(f"skip: line {t}") for l in log)
    assert sum(1 for l in log if l.startswith("ANN5320")) == 1 and sum(1 for l in log if l.startswith("ANN4300")) == 1


def test_apply_reports_missing_code():
    ann = _ann()
    log = apply(ann, [], {"ANN4300": ("misc", {})})
    assert log == ["no findings for ANN4300"] and _keys(ann) == ["source", "mRNA", "CDS", "mRNA", "CDS", "tRNA"]


def test_paired_mrna_with_spaces_in_ids(tmp_path):
    rows = ["NC_1\tsource\t1..100\torganism\tX"]
    for i in (1, 2):
        rows += [f"\tmRNA\t{i*10}..{i*10+5}\tlocus_tag\tL_{i}",
                 f"\t\t\tnote\tsubmitter_gene_id: Saccharomyces cerevisiae_NC_1_{i}, submitter_transcript_id: Saccharomyces cerevisiae_NC_1_{i}.1",
                 f"\tCDS\t{i*10}..{i*10+5}\tlocus_tag\tL_{i}",
                 f"\t\t\tnote\tsubmitter_gene_id: Saccharomyces cerevisiae_NC_1_{i}, submitter_transcript_id: Saccharomyces cerevisiae_NC_1_{i}.1"]
    p = tmp_path / "s.ann"
    p.write_text("\n".join(rows) + "\n")
    ann = AnnFile.read(p)
    cds2 = [f for f in ann.features() if f.key == "CDS"][1]
    m = paired_mrna(ann, cds2)
    assert m is not None and m.location == "20..25"
