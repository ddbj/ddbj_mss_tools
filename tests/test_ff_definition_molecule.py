"""Tests for mol_type-driven ff_definition molecule token."""

import pytest

from common.source_builder import _molecule_token, ff_definition, SequenceRoleEntry


@pytest.mark.parametrize("mol_type,expected", [
    ("genomic DNA", "DNA"),
    ("other DNA", "DNA"),
    ("unassigned DNA", "DNA"),
    ("mRNA", "mRNA"),
    ("tRNA", "tRNA"),
    ("rRNA", "rRNA"),
    ("genomic RNA", "RNA"),
    ("transcribed RNA", "RNA"),
    ("viral cRNA", "RNA"),
    ("protein", "DNA"),
    ("", "DNA"),
    (None, "DNA"),
])
def test_molecule_token(mol_type, expected):
    assert _molecule_token(mol_type) == expected


def test_ff_definition_unplaced_wgs_dna():
    out = ff_definition(None, "seq1", "Homo sapiens", "", "genomic DNA", is_wgs=True)
    assert out == "Homo sapiens DNA, seq1"


def test_ff_definition_unplaced_wgs_rna():
    out = ff_definition(None, "seq1", "Homo sapiens", "", "genomic RNA", is_wgs=True)
    assert out == "Homo sapiens RNA, seq1"


def test_ff_definition_unplaced_wgs_mrna():
    out = ff_definition(None, "seq1", "Homo sapiens", "", "mRNA", is_wgs=True)
    assert out == "Homo sapiens mRNA, seq1"


def test_ff_definition_chromosome_complete_rna():
    e = SequenceRoleEntry("seq1", "chromosome", "1", "complete", False)
    out = ff_definition(e, "seq1", "Homo sapiens", "strainX", "genomic RNA", is_wgs=False)
    assert out == "Homo sapiens strainX RNA, chromosome 1, complete sequence"


def test_ff_definition_empty_mol_type_defaults_dna():
    out = ff_definition(None, "seq1", "Homo sapiens", "", "", is_wgs=True)
    assert out == "Homo sapiens DNA, seq1"
