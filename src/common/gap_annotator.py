"""
Shared assembly_gap feature annotator for DDBJ MSS.
"""

from __future__ import annotations

import re
from typing import Optional

Row = list[str]

# Recommended (gap_type, estimated_length) per linkage_evidence
_RECOMMENDED: dict[str, tuple[str, str]] = {
    "paired-ends":        ("within scaffold", "known"),
    "proximity ligation": ("within scaffold", "unknown"),
    "align genus":        ("within scaffold", "unknown"),
}


class GapAnnotator:
    """Detect N-runs in a sequence and emit assembly_gap feature rows.

    Parameters
    ----------
    linkage_evidence:
        Evidence type for the gap (e.g. ``"paired-ends"``).  Must be a key
        of ``_RECOMMENDED`` unless *gap_type* and *estimated_length* are both
        supplied explicitly.
    min_gap_length:
        Minimum consecutive N count to annotate as a gap (default: 10).
    gap_type:
        Value for the ``gap_type`` qualifier.  ``None`` (default) applies
        the recommended value for *linkage_evidence*.
    estimated_length:
        Value for the ``estimated_length`` qualifier (``"known"``,
        ``"unknown"``, or ``None`` to use the recommendation).
    expand_meta_expression:
        When ``True`` and the resolved *estimated_length* is ``"known"``,
        the qualifier value is replaced with the actual integer gap length
        instead of the string ``"known"``.
    """

    def __init__(
        self,
        linkage_evidence: str,
        min_gap_length: int = 10,
        gap_type: Optional[str] = None,
        estimated_length: Optional[str] = None,
        expand_meta_expression: bool = False,
    ) -> None:
        self.linkage_evidence = linkage_evidence
        self.min_gap_length = min_gap_length
        self.expand_meta_expression = expand_meta_expression

        rec = _RECOMMENDED.get(linkage_evidence)

        if gap_type is None:
            if rec is None:
                raise ValueError(
                    f"No recommended gap_type for linkage_evidence={linkage_evidence!r}. "
                    "Specify gap_type explicitly."
                )
            self.gap_type = rec[0]
        else:
            self.gap_type = gap_type

        if estimated_length is None:
            if rec is None:
                raise ValueError(
                    f"No recommended estimated_length for linkage_evidence={linkage_evidence!r}. "
                    "Specify estimated_length explicitly."
                )
            self.estimated_length = rec[1]
        else:
            self.estimated_length = estimated_length

        self._pattern = re.compile(f"[Nn]{{{min_gap_length},}}")

    def annotate(self, seq: str, seq_name: str | None = None) -> list[Row]:
        """Return assembly_gap feature rows for all N-runs in *seq*.

        Parameters
        ----------
        seq:
            Nucleotide sequence string.
        seq_name:
            When provided, set the first column of the first feature row to
            this value (used in draft genomes where source is described in
            COMMON).
        """
        rows: list[Row] = []
        for match in self._pattern.finditer(seq):
            start = match.start() + 1  # 1-based
            end = match.end()
            if self.expand_meta_expression and self.estimated_length == "known":
                length_val = str(end - start + 1)
            else:
                length_val = self.estimated_length
            rows.append(["", "assembly_gap", f"{start}..{end}", "estimated_length", length_val])
            rows.append(["", "", "", "gap_type", self.gap_type])
            rows.append(["", "", "", "linkage_evidence", self.linkage_evidence])
        if rows and seq_name:
            rows[0][0] = seq_name
        return rows
