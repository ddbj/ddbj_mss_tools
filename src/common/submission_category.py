"""
Submission category manager for DDBJ MSS tools.

Loads category rules from submission_categories.json and provides helpers to
inject DATATYPE/DIVISION/KEYWORD defaults into a common_dict, and to expose
source-feature qualifiers that are automatically added for a given category
(e.g. environmental_sample for ENV/MAG/MAG-WGS).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

_DEFINITIONS_PATH = Path(__file__).parent / "data" / "submission_categories.json"
_definitions: dict | None = None


def _load_definitions() -> dict:
    global _definitions
    if _definitions is None:
        _definitions = json.loads(_DEFINITIONS_PATH.read_text(encoding="utf-8"))
    return _definitions


def _resolve(name: str, visiting: frozenset[str] = frozenset()) -> dict:
    """Return the fully-resolved category dict (inheritance applied, own fields win)."""
    definitions = _load_definitions()
    if name not in definitions:
        valid = sorted(definitions)
        raise ValueError(
            f"Unknown _submission_category: '{name}'. Valid values: {valid}"
        )
    if name in visiting:
        raise ValueError(f"Circular inheritance detected in submission category '{name}'")

    own = dict(definitions[name])
    parents = own.pop("inherits", [])

    merged: dict = {}
    for parent in parents:
        parent_resolved = _resolve(parent, visiting | {name})
        _merge_into(merged, parent_resolved)

    _merge_into(merged, own)
    return merged


def _merge_into(base: dict, override: dict) -> None:
    """Merge *override* into *base* in-place using category-specific merge rules."""
    scalar_keys = {
        "required_datatype", "required_division", "required_tagset_id", "source_identifier",
    }
    list_union_keys = {
        "required_st_comments", "required_moltype", "invalid_moltype",
        "required_dblinks", "required_source_qualifiers",
    }

    for key, value in override.items():
        if key in scalar_keys:
            base[key] = value
        elif key in list_union_keys:
            existing: list = list(base.get(key, []))
            for item in value:
                if item not in existing:
                    existing.append(item)
            base[key] = existing
        elif key == "required_keywords":
            # MERGE: extend parent keywords; dedup by exact inner-list match.
            existing = list(base.get(key, []))
            for group in value:
                if group not in existing:
                    existing.append(group)
            base[key] = existing
        elif key == "auto_source_qualifiers":
            base[key] = {**base.get(key, {}), **value}
        else:
            base[key] = value


@dataclass
class SubmissionCategoryRules:
    """Fully-resolved rules for one submission category."""

    datatype: str | None = None
    division: str | None = None
    required_keywords: list[list[str]] = field(default_factory=list)
    auto_source_qualifiers: dict = field(default_factory=dict)
    required_tagset_id: str | None = None
    required_st_comments: list[str] = field(default_factory=list)
    required_moltype: list[str] = field(default_factory=list)
    invalid_moltype: list[str] = field(default_factory=list)
    required_dblinks: list[str] = field(default_factory=list)
    required_source_qualifiers: list[str] = field(default_factory=list)
    source_identifier: str | None = None


def get_category_rules(name: str) -> SubmissionCategoryRules:
    """
    Return resolved rules for *name*. Raises ValueError for unknown categories.
    'NONE' (or empty string) is valid and returns an empty SubmissionCategoryRules.
    """
    if not name or name == "NONE":
        return SubmissionCategoryRules()
    resolved = _resolve(name)
    return SubmissionCategoryRules(
        datatype=resolved.get("required_datatype"),
        division=resolved.get("required_division"),
        required_keywords=resolved.get("required_keywords", []),
        auto_source_qualifiers=resolved.get("auto_source_qualifiers", {}),
        required_tagset_id=resolved.get("required_tagset_id"),
        required_st_comments=resolved.get("required_st_comments", []),
        required_moltype=resolved.get("required_moltype", []),
        invalid_moltype=resolved.get("invalid_moltype", []),
        required_dblinks=resolved.get("required_dblinks", []),
        required_source_qualifiers=resolved.get("required_source_qualifiers", []),
        source_identifier=resolved.get("source_identifier"),
    )


def validate_and_fill(common_dict: dict, category: str) -> None:
    """
    Check required fields; warn to stderr and fill with '' if missing.

    - required_source_qualifiers → common_dict["SOURCE"]
      (keys already in auto_source_qualifiers are skipped — they are added automatically)
    - required_dblinks           → common_dict["DBLINK"]
    - required_st_comments       → common_dict["ST_COMMENT"]
    """
    import sys
    rules = get_category_rules(category)
    auto_keys = set(rules.auto_source_qualifiers)

    user_required = [k for k in rules.required_source_qualifiers if k not in auto_keys]
    if user_required:
        source = common_dict.setdefault("SOURCE", {})
        for key in user_required:
            if not source.get(key):
                print(
                    f"[WARNING] [{category}] Required source qualifier '{key}' is missing."
                    " Adding empty value.",
                    file=sys.stderr,
                )
                source[key] = ""

    if rules.required_dblinks:
        dblink = common_dict.get("DBLINK", {})
        for key in rules.required_dblinks:
            # Accept both "sequence read archive" (original) and "sequence_read_archive" (pydantic serialized)
            alt_key = key.replace(" ", "_")
            if not dblink.get(key) and not dblink.get(alt_key):
                print(
                    f"[WARNING] [{category}] Required DBLINK '{key}' is missing."
                    " Adding empty value.",
                    file=sys.stderr,
                )
                common_dict.setdefault("DBLINK", {})[key] = ""

    if rules.required_st_comments:
        st = common_dict.get("ST_COMMENT")
        if st is None:
            print(
                f"[WARNING] [{category}] ST_COMMENT is missing."
                f" Required fields: {rules.required_st_comments}. Adding empty values.",
                file=sys.stderr,
            )
            common_dict["ST_COMMENT"] = {k: "" for k in rules.required_st_comments}
        else:
            for item in (st if isinstance(st, list) else [st]):
                for key in rules.required_st_comments:
                    if not item.get(key):
                        print(
                            f"[WARNING] [{category}] Required ST_COMMENT field '{key}' is missing."
                            " Adding empty value.",
                            file=sys.stderr,
                        )
                        item[key] = ""


def inject_defaults(common_dict: dict, category: str) -> None:
    """
    Inject DATATYPE / DIVISION / KEYWORD defaults into *common_dict* in-place.

    Only sets a field when it is absent or empty ({} / []).
    For KEYWORD, each required_keywords inner-list is checked: if none of its
    elements appear in the existing keyword list, the first element (the default)
    is appended.
    """
    rules = get_category_rules(category)

    if rules.datatype:
        current = common_dict.get("DATATYPE")
        if not current:
            common_dict["DATATYPE"] = {"type": rules.datatype}

    if rules.division:
        current = common_dict.get("DIVISION")
        if not current:
            common_dict["DIVISION"] = {"division": rules.division}

    if rules.required_keywords:
        existing_kw: list = list(common_dict.get("KEYWORD", {}).get("keyword", []))
        existing_set = set(existing_kw)
        to_add: list[str] = []
        for group in rules.required_keywords:
            if not any(k in existing_set for k in group):
                to_add.append(group[0])
        if to_add:
            common_dict["KEYWORD"] = {"keyword": existing_kw + to_add}
