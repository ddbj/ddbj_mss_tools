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


from common.source_builder import create_source_feature


def _ff_def_value(rows):
    for r in rows:
        if r[3] == "ff_definition":
            return r[4]
    return None


def test_create_source_meta_b2_rna():
    src = {"organism": "X", "mol_type": "genomic RNA"}
    rows = create_source_feature("WGS", None, None, None, src,
                                 source_modifier_key="strain",
                                 use_meta_expression=True)
    assert _ff_def_value(rows) == "@@[organism]@@ @@[strain]@@ RNA, @@[submitter_seqid]@@"


def test_create_source_meta_b2_default_dna():
    src = {"organism": "X"}
    rows = create_source_feature("WGS", None, None, None, src,
                                 source_modifier_key="strain",
                                 use_meta_expression=True)
    assert _ff_def_value(rows) == "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"


def test_create_source_b1_complete_rna():
    src = {"organism": "X", "mol_type": "genomic RNA"}
    rows = create_source_feature("GNM", "chr1", "complete", "linear", src,
                                 source_modifier_key="strain")
    assert _ff_def_value(rows) == "@@[organism]@@ @@[strain]@@ RNA, complete genome"


def test_create_source_b1_default_dna():
    src = {"organism": "X"}
    rows = create_source_feature("GNM", "chr1", "complete", "linear", src,
                                 source_modifier_key="strain")
    assert _ff_def_value(rows) == "@@[organism]@@ @@[strain]@@ DNA, complete genome"
