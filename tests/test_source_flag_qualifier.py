"""Tests for source flag-qualifier truthiness handling."""

from common.source_builder import _flag_is_set


def test_flag_is_set_false_strings():
    for v in ["false", "False", "FALSE", "no", "No", "NO"]:
        assert _flag_is_set(v) is False, f"{v!r} should be off"


def test_flag_is_set_false_bool():
    assert _flag_is_set(False) is False


def test_flag_is_set_true_bool():
    assert _flag_is_set(True) is True


def test_flag_is_set_truthy_strings():
    for v in ["yes", "Yes", "true", "True", "1", "x", "ATCC"]:
        assert _flag_is_set(v) is True, f"{v!r} should be on"


def test_flag_is_set_empty_string_is_on():
    # backward compatibility: JSON "environmental_sample": "" means on
    assert _flag_is_set("") is True


def test_flag_is_set_whitespace_false():
    assert _flag_is_set("  no  ") is False


from common.source_builder import _source_qualifier_rows


def test_rows_flag_on_emits_valueless_row():
    assert _source_qualifier_rows("environmental_sample", "yes") == [
        ["", "", "", "environmental_sample", ""]
    ]


def test_rows_flag_empty_string_on():
    assert _source_qualifier_rows("environmental_sample", "") == [
        ["", "", "", "environmental_sample", ""]
    ]


def test_rows_flag_off_emits_nothing():
    assert _source_qualifier_rows("environmental_sample", "no") == []
    assert _source_qualifier_rows("transgenic", False) == []


def test_rows_nonflag_keeps_value():
    # a non-flag qualifier whose value happens to be "No" is NOT treated as a flag
    assert _source_qualifier_rows("strain", "No") == [
        ["", "", "", "strain", "No"]
    ]


def test_rows_nonflag_list_multiple():
    assert _source_qualifier_rows("culture_collection", ["ATCC:1", "NBRC:2"]) == [
        ["", "", "", "culture_collection", "ATCC:1"],
        ["", "", "", "culture_collection", "NBRC:2"],
    ]


from common.models import CommonModel


def test_source_accepts_bool_value():
    data = {
        "DBLINK": {"project": "P", "sample": "S"},
        "SOURCE": {"organism": "E. coli", "environmental_sample": False},
    }
    m = CommonModel.model_validate(data)
    assert m.SOURCE["environmental_sample"] is False


def test_source_accepts_list_value_still():
    data = {
        "DBLINK": {"project": "P", "sample": "S"},
        "SOURCE": {"organism": "E. coli", "culture_collection": ["ATCC:1", "NBRC:2"]},
    }
    m = CommonModel.model_validate(data)
    assert m.SOURCE["culture_collection"] == ["ATCC:1", "NBRC:2"]


from common.source_builder import create_source_feature


def _quals(rows):
    """(key, value) pairs from 5-column source rows, skipping the feature line."""
    return [(r[3], r[4]) for r in rows if r[3]]


def test_create_source_flag_off_omitted_per_entry():
    src = {"organism": "E. coli", "environmental_sample": "no"}
    rows = create_source_feature("GNM", "chr1", "complete", "linear", src,
                                 source_modifier_key="strain")
    keys = [k for k, _ in _quals(rows)]
    assert "environmental_sample" not in keys
    assert "organism" in keys


def test_create_source_flag_on_valueless_per_entry():
    src = {"organism": "E. coli", "environmental_sample": "yes"}
    rows = create_source_feature("GNM", "chr1", "complete", "linear", src,
                                 source_modifier_key="strain")
    assert ("environmental_sample", "") in _quals(rows)


def test_create_source_nonflag_no_value_kept():
    # strain="No" must remain a normal valued qualifier, not be dropped
    src = {"organism": "E. coli", "strain": "No"}
    rows = create_source_feature("GNM", "chr1", "complete", "linear", src,
                                 source_modifier_key="strain")
    assert ("strain", "No") in _quals(rows)


def test_create_source_flag_on_meta_path():
    src = {"organism": "E. coli", "transgenic": True}
    rows = create_source_feature("WGS", None, None, None, src,
                                 source_modifier_key="strain",
                                 use_meta_expression=True)
    assert ("transgenic", "") in _quals(rows)


def test_create_source_flag_off_meta_path():
    src = {"organism": "E. coli", "transgenic": False}
    rows = create_source_feature("WGS", None, None, None, src,
                                 source_modifier_key="strain",
                                 use_meta_expression=True)
    keys = [k for k, _ in _quals(rows)]
    assert "transgenic" not in keys
