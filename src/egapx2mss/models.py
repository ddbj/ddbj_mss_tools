"""
Pydantic models for DDBJ MSS common metadata JSON.

The JSON file is loaded with load_common_json(), which also handles
trailing commas (JSON5-style relaxed syntax).
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, field_validator


class DblinkModel(BaseModel):
    project: str          # required
    sample: str           # required
    DRA: Optional[list[str]] = None


class SubmitterModel(BaseModel):
    ab_name: list[str]
    contact: Optional[str] = None
    email: Optional[str] = None
    phone: Optional[str] = None
    fax: Optional[str] = None
    institute: Optional[str] = None
    department: Optional[str] = None
    country: Optional[str] = None
    state: Optional[str] = None
    city: Optional[str] = None
    street: Optional[str] = None
    zip: Optional[str] = None

    @field_validator("ab_name")
    @classmethod
    def ab_name_not_empty(cls, v: list[str]) -> list[str]:
        if not v:
            raise ValueError("SUBMITTER.ab_name must contain at least one name")
        return v


class ReferenceModel(BaseModel):
    title: str
    ab_name: list[str]
    status: str
    year: int
    journal: Optional[str] = None
    volume: Optional[str] = None
    start_page: Optional[str] = None
    end_page: Optional[str] = None


_VALID_LINKAGE_EVIDENCE = {"paired-ends", "proximity ligation", "align genus"}


class AssemblyGapModel(BaseModel):
    linkage_evidence: str
    min_gap_length: int = 10

    @field_validator("linkage_evidence")
    @classmethod
    def check_linkage_evidence(cls, v: str) -> str:
        if v not in _VALID_LINKAGE_EVIDENCE:
            raise ValueError(
                f"linkage_evidence must be one of {sorted(_VALID_LINKAGE_EVIDENCE)}, got '{v}'"
            )
        return v


class CommonModel(BaseModel):
    DBLINK: DblinkModel
    SUBMITTER: Optional[SubmitterModel] = None
    REFERENCE: Optional[list[ReferenceModel]] = None
    DATE: Optional[dict[str, str]] = None
    SOURCE: Optional[dict[str, str]] = None
    DATATYPE: Optional[dict[str, str]] = None
    KEYWORD: Optional[dict[str, str | list[str]]] = None
    ASSEMBLY_GAP: Optional[AssemblyGapModel] = None
    SOURCE_IDENTIFIER: Optional[str] = None


def _strip_trailing_commas(text: str) -> str:
    """Remove trailing commas before closing braces/brackets (JSON5 relaxation)."""
    return re.sub(r",\s*([}\]])", r"\1", text)


def load_common_json(path: str) -> CommonModel:
    """
    Load and validate a common metadata JSON file.

    Trailing commas before ``}`` or ``]`` are accepted (JSON5-style).
    Raises ``pydantic.ValidationError`` if required fields are missing.
    """
    raw = Path(path).read_text(encoding="utf-8")
    data = json.loads(_strip_trailing_commas(raw))
    return CommonModel.model_validate(data)
