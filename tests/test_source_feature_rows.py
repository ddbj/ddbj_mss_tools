"""Per-entry source feature emission must apply the same flag/bool/list handling
as the COMMON path (via _source_qualifier_rows), not dump raw common.SOURCE
values into rows. Regression for: a boolean flag qualifier (e.g.
"environmental_sample": true) crashing mss_builder with
TypeError: sequence item ...: expected str instance, bool found.
"""

from common.models import CommonModel
from common.source_builder import source_feature_rows, SequenceRoleEntry
from mss_builder.ann_writer import write_mss_ann


# ── unit: source_feature_rows helper ───────────────────────────────────────
def test_source_feature_rows_bool_flag_on_valueless():
    rows = source_feature_rows("chr1", "1..100", {
        "organism": "E. coli", "mol_type": "genomic DNA", "environmental_sample": True})
    assert rows[0][:3] == ["chr1", "source", "1..100"]          # first row is the header
    assert ["", "", "", "environmental_sample", ""] in rows      # valueless flag row
    assert all(isinstance(c, str) for r in rows for c in r)      # no bool leaked


def test_source_feature_rows_bool_flag_off_omitted():
    rows = source_feature_rows("chr1", "1..100", {
        "organism": "E. coli", "environmental_sample": False})
    assert not any(r[3] == "environmental_sample" for r in rows)
    assert all(isinstance(c, str) for r in rows for c in r)


def test_source_feature_rows_list_value_expands():
    rows = source_feature_rows("chr1", "1..100", {
        "organism": "E. coli", "culture_collection": ["ATCC:1", "NBRC:2"]})
    assert [r[4] for r in rows if r[3] == "culture_collection"] == ["ATCC:1", "NBRC:2"]


def test_source_feature_rows_first_qualifier_is_header():
    rows = source_feature_rows("c", "1..10", {"organism": "X", "strain": "s1"})
    assert rows[0] == ["c", "source", "1..10", "organism", "X"]
    assert rows[1] == ["", "", "", "strain", "s1"]


# ── integration: write_mss_ann non-WGS per-entry path (the crash scenario) ──
def _common_with(env):
    return CommonModel.model_validate({
        "DBLINK": {"project": "P", "sample": "S"},
        "SOURCE": {"organism": "Bos taurus", "mol_type": "genomic DNA",
                   "strain": "Example-1", "environmental_sample": env},
    })


def _run(tmp_path, env):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">chr1\nACGTACGTAC\n")
    ann = tmp_path / "out.ann"
    roles = {"chr1": SequenceRoleEntry("chr1", "chromosome", "1", "complete", False)}
    write_mss_ann(str(fasta), str(ann), common=_common_with(env), sequence_roles=roles)
    return ann.read_text()


def test_write_mss_ann_bool_flag_true_per_entry_no_crash(tmp_path):
    text = _run(tmp_path, True)                    # must not raise TypeError
    assert "\t\t\tenvironmental_sample\t\n" in text  # valueless flag row emitted


def test_write_mss_ann_bool_flag_false_per_entry_omitted(tmp_path):
    text = _run(tmp_path, False)                    # must not raise TypeError
    assert "environmental_sample" not in text        # flag row omitted when off
