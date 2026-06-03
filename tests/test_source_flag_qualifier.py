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
