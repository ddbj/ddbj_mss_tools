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


from common.source_builder import create_source_feature, source_qualifier, _organelle_code


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


# ── organelle code helper ──────────────────────────────────────────────
@pytest.mark.parametrize("seq_name,expected", [
    ("mitochondrion", "mitochondrial"),
    ("mitochondrion:kinetoplast", "kinetoplast"),
    ("hydrogenosome", "hydrogenosomal"),
    ("nucleomorph", "nucleomorph"),
    ("plastid", "plastid"),
    ("plastid:chloroplast", "chloroplast"),
    ("plastid:apicoplast", "apicoplast"),
    ("plastid:chromoplast", "chromoplast"),
    ("plastid:cyanelle", "cyanelle"),
    ("plastid:leucoplast", "leucoplast"),
    ("plastid:proplastid", "proplastid"),
    ("macronuclear", "macronuclear"),
    ("mitochondrial", "mitochondrial"),   # not in table -> passthrough
    ("anything else", "anything else"),
])
def test_organelle_code(seq_name, expected):
    assert _organelle_code(seq_name) == expected


# ── chromosome: count-dependent completeness ───────────────────────────
def test_ff_chromosome_single_complete_genome():
    e = SequenceRoleEntry("seq1", "chromosome", "Y", "complete", False)
    out = ff_definition(e, "seq1", "Pan troglodytes", "", "genomic DNA",
                        is_wgs=False, chromosome_count=1)
    assert out == "Pan troglodytes DNA, chromosome Y, complete genome"


def test_ff_chromosome_multi_complete_sequence():
    e = SequenceRoleEntry("seq1", "chromosome", "1", "complete", False)
    out = ff_definition(e, "seq1", "Homo sapiens", "", "genomic DNA",
                        is_wgs=False, chromosome_count=2)
    assert out == "Homo sapiens DNA, chromosome 1, complete sequence"


def test_ff_chromosome_default_count_complete_sequence():
    e = SequenceRoleEntry("seq1", "chromosome", "1", "complete", False)
    out = ff_definition(e, "seq1", "Homo sapiens", "", "genomic DNA", is_wgs=False)
    assert out == "Homo sapiens DNA, chromosome 1, complete sequence"


def test_ff_chromosome_partial_unlocalized_unchanged():
    e = SequenceRoleEntry("seq1", "chromosome", "1", "partial", False)
    out = ff_definition(e, "seq1", "Homo sapiens", "", "genomic DNA",
                        is_wgs=False, chromosome_count=1)
    assert out == "Homo sapiens DNA, chromosome 1, unlocalized sequence seq1"


# ── organelle: doc-compliant form ──────────────────────────────────────
def test_ff_organelle_mitochondrion_complete():
    e = SequenceRoleEntry("mt", "organelle", "mitochondrion", "complete", True)
    out = ff_definition(e, "mt", "Homo sapiens", "isolate TSX", "genomic DNA")
    assert out == "Homo sapiens isolate TSX mitochondrial DNA, complete genome"


def test_ff_organelle_partial_genome():
    e = SequenceRoleEntry("mt", "organelle", "mitochondrion", "partial", True)
    out = ff_definition(e, "mt", "Homo sapiens", "", "genomic DNA")
    assert out == "Homo sapiens mitochondrial DNA, partial genome"


def test_ff_organelle_plastid_chloroplast():
    e = SequenceRoleEntry("cp", "organelle", "plastid:chloroplast", "complete", True)
    out = ff_definition(e, "cp", "Zea mays", "", "genomic DNA")
    assert out == "Zea mays chloroplast DNA, complete genome"


def test_ff_organelle_passthrough_name():
    e = SequenceRoleEntry("mt", "organelle", "mitochondrial", "complete", True)
    out = ff_definition(e, "mt", "Homo sapiens", "", "genomic DNA")
    assert out == "Homo sapiens mitochondrial DNA, complete genome"


def test_ff_organelle_rna_token_position():
    e = SequenceRoleEntry("mt", "organelle", "mitochondrion", "complete", True)
    out = ff_definition(e, "mt", "Homo sapiens", "", "genomic RNA")
    assert out == "Homo sapiens mitochondrial RNA, complete genome"


# ── plasmid: new type ──────────────────────────────────────────────────
def test_ff_plasmid_complete():
    e = SequenceRoleEntry("p1", "plasmid", "pLG1", "complete", True)
    out = ff_definition(e, "p1", "Lactobacillus gasseri", "SG162", "genomic DNA")
    assert out == "Lactobacillus gasseri SG162 plasmid pLG1 DNA, complete sequence"


def test_ff_plasmid_partial():
    e = SequenceRoleEntry("p1", "plasmid", "pLG1", "partial", True)
    out = ff_definition(e, "p1", "Lactobacillus gasseri", "SG162", "genomic DNA")
    assert out == "Lactobacillus gasseri SG162 plasmid pLG1 DNA, partial sequence"


# ── source_qualifier: plasmid emitted, organelle stays RAW ─────────────
def test_source_qualifier_plasmid():
    e = SequenceRoleEntry("p1", "plasmid", "pLG1", "complete", True)
    assert source_qualifier(e, "p1") == {"plasmid": "pLG1"}


def test_source_qualifier_organelle_raw():
    e = SequenceRoleEntry("cp", "organelle", "plastid:chloroplast", "complete", True)
    assert source_qualifier(e, "cp") == {"organelle": "plastid:chloroplast"}
