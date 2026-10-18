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


# ---- review fix pass: identity by row object, feature-key guards, multi-id partners ----

def _write(tmp_path, rows):
    p = tmp_path / "t.ann"
    p.write_text("\n".join(rows) + "\n")
    return p


def _pair_rows(i, loc, tid, gene="G"):
    return [f"\tmRNA\t{loc}\tlocus_tag\tL_{i}",
            f"\t\t\tnote\tsubmitter_gene_id: {gene}, submitter_transcript_id: {tid}",
            f"\tCDS\t{loc}\tlocus_tag\tL_{i}",
            f"\t\t\tproduct\thypothetical protein",
            f"\t\t\tnote\tsubmitter_gene_id: {gene}, submitter_transcript_id: {tid}"]


def test_same_line_duplicate_findings_do_not_hit_next_feature(tmp_path):
    rows = ["NC_1\tsource\t1..1000\torganism\tX",
            "\ttRNA\t600..671\tlocus_tag\tL_1", "\t\t\tproduct\ttRNA-Pro",
            "\ttRNA\t700..771\tlocus_tag\tL_2", "\t\t\tproduct\ttRNA-Gly"]
    ann = AnnFile.read(_write(tmp_path, rows))
    msg = "... (Found: tRNA 600..671 in CDS 500..800)"
    findings = [Finding("ANN5320", "ERR", "NC_1", "tRNA", 2, msg), Finding("ANN5320", "ERR", "NC_1", "tRNA", 2, msg)]
    log = apply(ann, findings, {"ANN5320": ("overlap-drop", {"keep": "cds"})})
    assert _keys(ann) == ["source", "tRNA"] and next(f for f in ann.features() if f.key == "tRNA").location == "700..771"
    assert sum(1 for l in log if l.startswith("skip")) == 1


def test_nonredundant_layout_two_cds_findings_in_one_gene(tmp_path):
    rows = ["NC_1\tsource\t1..1000\torganism\tX",
            "\tmRNA\t10..100\tlocus_tag\tL_1", "\t\t\tnote\tsubmitter_gene_id: G, submitter_transcript_id: T1",
            "\tmRNA\t10..200\tlocus_tag\tL_1", "\t\t\tnote\tsubmitter_gene_id: G, submitter_transcript_id: T2",
            "\tCDS\t10..100\tlocus_tag\tL_1", "\t\t\tproduct\tp", "\t\t\tnote\tsubmitter_gene_id: G, submitter_transcript_id: T1",
            "\tCDS\t10..200\tlocus_tag\tL_1", "\t\t\tproduct\tp", "\t\t\tnote\tsubmitter_gene_id: G, submitter_transcript_id: T2"]
    ann = AnnFile.read(_write(tmp_path, rows))
    findings = [Finding("ANN4300", "WAR", "NC_1", "CDS", 6, "Complete CDS must be at least 90 bases"),
                Finding("ANN4300", "WAR", "NC_1", "CDS", 9, "Complete CDS must be at least 90 bases")]
    apply(ann, findings, {"ANN4300": ("misc", {})})
    assert _keys(ann) == ["source", "misc_feature", "misc_feature"]


def test_identical_cds_isoforms_target_the_right_one(tmp_path):
    rows = ["NC_1\tsource\t1..1000\torganism\tX"] + _pair_rows(1, "10..100", "T1") + _pair_rows(1, "10..100", "T2")
    ann = AnnFile.read(_write(tmp_path, rows))
    cds_t2 = [f for f in ann.features() if f.key == "CDS"][1]
    apply(ann, [Finding("ANN4300", "WAR", "NC_1", "CDS", cds_t2.first + 1, "short")], {"ANN4300": ("misc", {})})
    feats = [f for f in ann.features() if f.entry != "COMMON"]
    assert [f.key for f in feats] == ["source", "mRNA", "CDS", "misc_feature"]
    assert "T1" in feats[1].get(ann.rows, "note") and "T2" in feats[3].qualifiers(ann.rows)[-1][1]


def test_shared_cds_removes_every_listed_mrna(tmp_path):
    rows = ["NC_1\tsource\t1..1000\torganism\tX",
            "\tmRNA\t10..100\tlocus_tag\tL_1", "\t\t\tnote\tsubmitter_gene_id: G, submitter_transcript_id: T1",
            "\tmRNA\t5..100\tlocus_tag\tL_1", "\t\t\tnote\tsubmitter_gene_id: G, submitter_transcript_id: T2",
            "\tCDS\t10..100\tlocus_tag\tL_1", "\t\t\tproduct\tp", "\t\t\tnote\tsubmitter_gene_id: G, submitter_transcript_id: T1, T2"]
    ann = AnnFile.read(_write(tmp_path, rows))
    cds = next(f for f in ann.features() if f.key == "CDS")
    op_misc(ann, cds, reason="x")
    assert _keys(ann) == ["source", "misc_feature"]


def test_artificial_location_skips_non_cds_features(tmp_path):
    rows = ["NC_1\tsource\t1..1000\torganism\tX", "\tintron\t50..60\t\t"]
    ann = AnnFile.read(_write(tmp_path, rows))
    log = apply(ann, [Finding("AXS6810", "WAR", "NC_1", "intron", 2, "Non-canonical splice sites (Found: GT-GT at 50..60)")],
                {"AXS6810": ("artificial-location", {"value": "low-quality sequence region"})})
    assert not any(k == "artificial_location" for f in ann.features() for k, _ in f.qualifiers(ann.rows))
    assert any("skip" in l and "intron" in l for l in log)


def test_drop_and_misc_refuse_source(tmp_path):
    rows = ["NC_1\tsource\t1..1000\torganism\tX", "\t\t\tstrain\tS288C"]
    ann = AnnFile.read(_write(tmp_path, rows))
    log = apply(ann, [Finding("ANN1110", "WAR", "NC_1", "source", 1, "strain matches an institution code")],
                {"ANN1110": ("drop", {})})
    assert _keys(ann) == ["source"] and any("skip" in l for l in log)


def test_unresolvable_finding_aborts_without_edits(tmp_path):
    ann = _ann()
    cds = _cds2(ann)
    before = [r.cols[:] for r in ann.rows]
    findings = [Finding("ANN4300", "WAR", "NC_1", "CDS", cds.first + 1, "short"),
                Finding("ANN4300", "WAR", "NC_1", "CDS", cds.first + 2, "short")]      # points at a qualifier row
    with pytest.raises(ValueError, match="does not match"):
        apply(ann, findings, {"ANN4300": ("misc", {})})
    assert [r.cols for r in ann.rows] == before


def test_location_in_message_must_match_the_row(tmp_path):
    ann = _ann()
    cds = _cds2(ann)
    with pytest.raises(ValueError, match="does not match"):
        apply(ann, [Finding("AXS6810", "WAR", "NC_1", "CDS", cds.first + 1, "(Found: GT-AT in CDS at join(1..2,3..4))")],
              {"AXS6810": ("artificial-location", {"value": "low-quality sequence region"})})


def test_keep_other_with_missing_cds_is_skipped_not_raised(tmp_path):
    ann = _ann()
    t = _trna(ann).first + 1
    log = apply(ann, [Finding("ANN5320", "ERR", "NC_1", "tRNA", t, "(Found: tRNA 600..671 in CDS 999..1000)")],
                {"ANN5320": ("overlap-drop", {"keep": "other"})})
    assert any(l.startswith("skip") for l in log) and _keys(ann) == ["source", "mRNA", "CDS", "mRNA", "CDS", "tRNA"]


def test_defaults_include_short_intron_code():
    assert DEFAULTS["ANN6840"] == ("artificial-location", {"value": "low-quality sequence region"})
