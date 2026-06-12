"""Tests for mol_type-driven ff_definition molecule token."""

import pytest

from common.source_builder import _molecule_token


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
