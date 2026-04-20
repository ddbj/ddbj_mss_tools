"""
Pydantic models for DDBJ MSS common metadata JSON.

Re-exported from common.models for backwards compatibility.
New code should import directly from common.models.
"""

from common.models import (  # noqa: F401
    AssemblyGapModel,
    CommonModel,
    DblinkModel,
    ReferenceModel,
    SubmitterModel,
    _strip_trailing_commas,
    load_common_json,
)
