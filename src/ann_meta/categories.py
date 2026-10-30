"""Submission-category sets shared by ``validate``, ``apply`` and the CLI.

``category`` here is the mss_tools name passed to ``ann_meta apply
--category`` (a key of ``submission_categories.json``, or ``NONE`` / ``GENE``).

* ``WGS_LIKE`` -- contig sets: when every entry is unplaced without per-entry
  qualifiers, one COMMON source at ``1..E`` with ``submitter_seqid
  @@[entry]@@`` (mss_builder's WGS, tsa2mss's TSA); the ff_definition names
  the contig by ``@@[submitter_seqid]@@``.
* ``ROLE_CATEGORIES`` -- entries may carry a role (chromosome / plasmid /
  organelle / segment); every other category takes unplaced entries only.
* ``REGENERATE_FF`` -- ff_definition is rebuilt from the role
  (``source_builder.ff_definition``); elsewhere it is written as given.
* ``PASS_THROUGH`` -- gene-type / unknown submissions: roles, names and
  topology are written as given and the role rules are not checked (a
  plasmid-only or virus submission is circular and complete).
"""
from __future__ import annotations

WGS_LIKE = frozenset({"WGS", "MAG-WGS", "ENV", "TSA", "TSA-ENV"})
PASS_THROUGH = frozenset({"GENE", "NONE"})
ROLE_CATEGORIES = frozenset({"GNM", "MAG", "TPA-GNM"}) | PASS_THROUGH
REGENERATE_FF = frozenset({"GNM", "WGS", "MAG", "MAG-WGS", "ENV", "TSA", "TSA-ENV"})
