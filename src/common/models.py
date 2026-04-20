"""
Pydantic models for DDBJ MSS common metadata JSON.

Shared between egapx2mss and mss_builder.
The JSON file is loaded with load_common_json(), which also handles
trailing commas (JSON5-style relaxed syntax).
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, field_validator, model_validator


class DblinkModel(BaseModel):
    project: str          # required
    sample: str           # required
    DRA: Optional[list[str]] = None


class SubmitterModel(BaseModel):
    ab_name: list[str] = []
    consrtm: Optional[str] = None
    contact: Optional[str] = None
    email: Optional[str] = None
    phone: Optional[str] = None
    fax: Optional[str] = None
    url: Optional[str] = None
    institute: Optional[str] = None
    department: Optional[str] = None
    country: Optional[str] = None
    state: Optional[str] = None
    city: Optional[str] = None
    street: Optional[str] = None
    zip: Optional[str] = None

    @model_validator(mode="after")
    def require_author_or_consrtm(self) -> "SubmitterModel":
        if not self.ab_name and not self.consrtm:
            raise ValueError("SUBMITTER must have at least one ab_name or a consrtm")
        return self


_VALID_STATUSES = {"Unpublished", "In press", "Published"}


class ReferenceModel(BaseModel):
    title: str
    ab_name: list[str] = []
    consrtm: Optional[str] = None
    status: str
    year: Optional[int] = None
    journal: Optional[str] = None
    volume: Optional[str] = None
    start_page: Optional[str] = None
    end_page: Optional[str] = None

    @model_validator(mode="after")
    def validate_status_fields(self) -> "ReferenceModel":
        if not self.ab_name and not self.consrtm:
            raise ValueError("REFERENCE must have at least one ab_name or a consrtm")
        status = self.status
        if status not in _VALID_STATUSES:
            raise ValueError(
                f"REFERENCE.status must be one of {sorted(_VALID_STATUSES)}, got '{status}'"
            )
        if status == "In press":
            if not self.journal:
                raise ValueError("REFERENCE.journal is required when status='In press'")
            if self.year is None:
                raise ValueError("REFERENCE.year is required when status='In press'")
        elif status == "Published":
            missing = [f for f in ("journal", "volume", "start_page") if not getattr(self, f)]
            if missing:
                raise ValueError(
                    f"REFERENCE fields required when status='Published': {', '.join(missing)}"
                )
            if self.year is None:
                raise ValueError("REFERENCE.year is required when status='Published'")
        return self


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
    SOURCE_MODIFIER: Optional[str] = None


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
