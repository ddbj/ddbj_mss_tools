from pathlib import Path

from ann_fix.report import Finding, group_by_code, parse_details

FX = Path(__file__).parent / "ann_fix_fixtures" / "details.txt"


def test_counts_per_code():
    groups = group_by_code(parse_details(FX.read_text()))
    assert {c: len(v) for c, v in groups.items()} == {
        "ANN5320": 2, "ANN4300": 2, "AXS6810": 2, "ANN2520": 1, "ANN2542": 1, "ANN0470": 1, "ANN2670": 1}
    assert list(groups)[:2] == ["ANN5320", "ANN4300"]        # first-appearance order


def test_fields_and_levels():
    f = group_by_code(parse_details(FX.read_text()))
    e = f["ANN5320"][0]
    assert e == Finding("ANN5320", "ERR", "NC_001142.9", "tRNA", 22138,
                        e.message) and e.message.startswith("The tRNA feature must not")
    assert f["ANN2520"][0].line == 1405 and f["ANN2520"][0].feature == "CDS"
    assert f["ANN2542"][0].line is None and f["ANN2542"][0].entry == "ALL"
    assert f["ANN0470"][0].line is None
    assert f["ANN2670"][0].level == "INF" and f["ANN2670"][0].line == 30


def test_headers_and_blank_lines_ignored():
    assert parse_details("=== Validation Details ===\nVersion: 1.0.70\n\n") == []
