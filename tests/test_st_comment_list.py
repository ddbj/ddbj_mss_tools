"""ST_COMMENT tag values may be a list, joined with '; ' into a single row.
Strings (including ones that already contain '; ') pass through unchanged.
Other features keep list == one-row-per-element behavior.
"""

from common.common_builder import create_feature


def _val(rows, key):
    """Values emitted for a given qualifier key (one entry per output row)."""
    return [r[4] for r in rows if r[3] == key]


def test_st_comment_list_joined_with_semicolon():
    rows = create_feature("ST_COMMENT", {
        "tagset_id": "Genome-Assembly-Data",
        "Assembly Method": ["Skesa v. 1.0", "Hifiasm", "yahs v. 2.1"],
        "Sequencing Technology": ["Illumina NanoSeq", "PacBio Revio"],
    })
    # each list collapses to a SINGLE row joined with '; '
    assert _val(rows, "Assembly Method") == ["Skesa v. 1.0; Hifiasm; yahs v. 2.1"]
    assert _val(rows, "Sequencing Technology") == ["Illumina NanoSeq; PacBio Revio"]


def test_st_comment_string_with_semicolon_passthrough():
    rows = create_feature("ST_COMMENT", {
        "tagset_id": "Genome-Assembly-Data",
        "Sequencing Technology": "Illumina NanoSeq; PacBio Revio",
    })
    assert _val(rows, "Sequencing Technology") == ["Illumina NanoSeq; PacBio Revio"]


def test_st_comment_single_string_unchanged():
    rows = create_feature("ST_COMMENT", {
        "tagset_id": "X", "Assembly Method": "SPAdes v. 3.15.5"})
    assert _val(rows, "Assembly Method") == ["SPAdes v. 3.15.5"]


def test_st_comment_single_element_list():
    rows = create_feature("ST_COMMENT", {"tagset_id": "X", "Assembly Method": ["Hifiasm"]})
    assert _val(rows, "Assembly Method") == ["Hifiasm"]


def test_non_st_comment_list_still_multiple_rows():
    # KEYWORD (and REFERENCE/SOURCE etc.) keep list == one row per element
    rows = create_feature("KEYWORD", {"keyword": ["WGS", "STANDARD_DRAFT"]})
    assert _val(rows, "keyword") == ["WGS", "STANDARD_DRAFT"]
