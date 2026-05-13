"""
Shared utilities for writing DDBJ MSS annotation files from JSON metadata.

Used by both egapx2mss and wgs_maker.
"""

import json

from common.source_builder import create_source_feature

# Keys in the common JSON that are tool-specific configuration, not DDBJ MSS features
_NON_COMMON_KEYS = frozenset({"SOURCE", "SOURCE_IDENTIFIER", "INFRASPECIFIC_NAME_MODIFIER", "ASSEMBLY_GAP"})


def create_qualifier(qualifier_key: str, value: str | list) -> list[list[str]]:
    ret = []
    if isinstance(value, list):
        for v in value:
            ret.append(["", "", "", qualifier_key, str(v)])
    else:
        ret.append(["", "", "", qualifier_key, str(value)])
    return ret


def sort_st_comment(qualifier_keys: list[str]) -> list[str]:
    return [k for k in qualifier_keys if k == "tagset_id"] + sorted(
        [k for k in qualifier_keys if k != "tagset_id"]
    )


def create_feature(feature_name: str, feature_values: dict | list) -> list[list[str]]:
    ret = []
    if isinstance(feature_values, list):  # REFERENCE and COMMENT are arrays
        for v in feature_values:
            ret.extend(create_feature(feature_name, v))
    elif feature_name == "ST_COMMENT":  # qualifiers must be sorted alphabetically
        qualifier_keys = sort_st_comment(feature_values.keys())
        for qualifier_key in qualifier_keys:
            ret.extend(create_qualifier(qualifier_key, feature_values[qualifier_key]))
    else:
        for qualifier_key, value in feature_values.items():
            ret.extend(create_qualifier(qualifier_key, value))
    if ret:
        ret[0][1] = feature_name
    return ret


def _build_common_source(common_json: dict) -> list[list[str]]:
    """Build source feature rows with meta-notation for inclusion in the COMMON block.

    Delegates to :func:`source_builder.create_source_feature` with
    ``use_meta_expression=True``.  The source qualifiers, submission category, and
    ff_definition modifier key are read from the ``SOURCE``, ``_trad_submission_category``,
    and ``SOURCE_IDENTIFIER`` keys of *common_json*.
    """
    source_data: dict = common_json.get("SOURCE", {})
    infraspecific_name_modifier_key: str = (
        common_json.get("SOURCE_IDENTIFIER") or common_json.get("INFRASPECIFIC_NAME_MODIFIER", "")
    )
    category: str = common_json.get("_trad_submission_category", "")
    return create_source_feature(
        category,
        None, None, None,
        source_data,
        source_modifier_key=infraspecific_name_modifier_key,
        use_meta_expression=True,
    )


def create_common(common_json: dict, include_source: bool = False) -> list[list[str]]:
    """Convert a common metadata dict to a list of MSS annotation rows (5-element lists).

    Keys starting with '_' and tool-specific keys (SOURCE, SOURCE_IDENTIFIER, ASSEMBLY_GAP)
    are skipped. The DBLINK qualifier 'sample' is output as 'biosample'.

    If *include_source* is True a source feature is appended to the COMMON block using
    ``1..E`` location and ``@@[...]@@`` meta-notation for per-entry fields.  This is
    intended for WGS submissions where a single COMMON source feature applies to all
    contigs.  The source qualifiers and ff_definition modifier are taken from the
    ``SOURCE`` and ``SOURCE_IDENTIFIER`` keys of *common_json*.
    """
    ret = []
    for feature_name, feature_values in common_json.items():
        if feature_name.startswith("_") or feature_name in _NON_COMMON_KEYS:
            continue
        if feature_name == "DBLINK" and isinstance(feature_values, dict):
            _dblink_key_map = {"sample": "biosample", "sequence_read_archive": "sequence read archive"}
            feature_values = {
                _dblink_key_map.get(k, k): v
                for k, v in feature_values.items()
            }
        ret.extend(create_feature(feature_name, feature_values))
    if include_source:
        ret.extend(_build_common_source(common_json))
    if ret:
        ret[0][0] = "COMMON"
    return ret


if __name__ == "__main__":
    import sys

    input_json = sys.argv[1]
    common_json = json.load(open(input_json))
    common = create_common(common_json)
    for row in common:
        print("\t".join(map(str, row)))
