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


# ── source_qualifier: plasmid emitted, organelle stays RAW ─────────────
def test_source_qualifier_plasmid():
    e = SequenceRoleEntry("p1", "plasmid", "pLG1", "complete", True)
    assert source_qualifier(e, "p1") == {"plasmid": "pLG1"}


def test_source_qualifier_organelle_raw():
    e = SequenceRoleEntry("cp", "organelle", "plastid:chloroplast", "complete", True)
    assert source_qualifier(e, "cp") == {"organelle": "plastid:chloroplast"}


# ── source_qualifier: segment emitted only when multi-segment ──────────────
def test_source_qualifier_segment_single_omitted():
    e = SequenceRoleEntry("seg1", "segment", "", "complete", False)
    assert source_qualifier(e, "seg1", segment_count=1) == {}


def test_source_qualifier_segment_single_name_ignored():
    e = SequenceRoleEntry("seg1", "segment", "4", "complete", False)
    assert source_qualifier(e, "seg1", segment_count=1) == {}


def test_source_qualifier_segment_multi():
    e = SequenceRoleEntry("seg4", "segment", "4", "complete", False)
    assert source_qualifier(e, "seg4", segment_count=8) == {"segment": "4"}


def test_source_qualifier_segment_multi_empty_name_omitted():
    e = SequenceRoleEntry("seg4", "segment", "", "complete", False)
    assert source_qualifier(e, "seg4", segment_count=8) == {}


def test_source_qualifier_segment_default_count_omitted():
    e = SequenceRoleEntry("seg1", "segment", "4", "complete", False)
    assert source_qualifier(e, "seg1") == {}


# ── ff_definition: meta-notation (new signature) ───────────────────────
def _entry(type_, seq_name="", status="complete"):
    return SequenceRoleEntry("sid", type_, seq_name, status, False)


@pytest.mark.parametrize("entry,source_identifier,mol_type,is_wgs,chrom,seg,expected", [
    # unplaced / None
    (None, None, "genomic DNA", True, 0, 0,
     "@@[organism]@@ DNA, @@[submitter_seqid]@@"),
    (None, None, "genomic DNA", False, 0, 0,
     "@@[organism]@@ DNA, unplaced sequence @@[entry]@@"),
    (None, "cultivar", "genomic RNA", False, 0, 0,
     "@@[organism]@@ @@[cultivar]@@ RNA, unplaced sequence @@[entry]@@"),
    (_entry("unplaced"), "strain", "genomic DNA", True, 0, 0,
     "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"),
    # chromosome — single (count<=1): no number
    (_entry("chromosome", "1", "complete"), "strain", "genomic DNA", False, 1, 0,
     "@@[organism]@@ @@[strain]@@ DNA, chromosome, complete genome"),
    (_entry("chromosome", "1", "partial"), "strain", "genomic DNA", False, 1, 0,
     "@@[organism]@@ @@[strain]@@ DNA, chromosome"),
    (_entry("chromosome", "", "complete"), None, "genomic DNA", False, 1, 0,
     "@@[organism]@@ DNA, chromosome, complete genome"),
    # chromosome — multiple (count>=2): @@[chromosome]@@
    (_entry("chromosome", "1", "complete"), "strain", "genomic DNA", False, 2, 0,
     "@@[organism]@@ @@[strain]@@ DNA, chromosome @@[chromosome]@@, complete sequence"),
    (_entry("chromosome", "1", "partial"), "strain", "genomic DNA", False, 2, 0,
     "@@[organism]@@ @@[strain]@@ DNA, chromosome @@[chromosome]@@"),
    # organelle — prefix meta, adjective concrete
    (_entry("organelle", "mitochondrion", "complete"), "", "genomic DNA", False, 0, 0,
     "@@[organism]@@ mitochondrial DNA, complete genome"),
    (_entry("organelle", "mitochondrion", "partial"), None, "genomic DNA", False, 0, 0,
     "@@[organism]@@ mitochondrial DNA, partial genome"),
    (_entry("organelle", "plastid:chloroplast", "complete"), "isolate", "genomic DNA", False, 0, 0,
     "@@[organism]@@ @@[isolate]@@ chloroplast DNA, complete genome"),
    (_entry("organelle", "mitochondrion", "complete"), None, "genomic RNA", False, 0, 0,
     "@@[organism]@@ mitochondrial RNA, complete genome"),
    # plasmid
    (_entry("plasmid", "pLG1", "complete"), "strain", "genomic DNA", False, 0, 0,
     "@@[organism]@@ @@[strain]@@ plasmid @@[plasmid]@@ DNA, complete sequence"),
    (_entry("plasmid", "pLG1", "partial"), "strain", "genomic DNA", False, 0, 0,
     "@@[organism]@@ @@[strain]@@ plasmid @@[plasmid]@@ DNA, partial sequence"),
    # segment — single (count<=1): no 'segment' word
    (_entry("segment", "", "complete"), "strain", "viral cRNA", False, 0, 1,
     "@@[organism]@@ @@[strain]@@ RNA, complete genome"),
    (_entry("segment", "", "partial"), "strain", "viral cRNA", False, 0, 1,
     "@@[organism]@@ @@[strain]@@ RNA, partial genome"),
    # segment — multiple (count>=2): @@[segment]@@
    (_entry("segment", "4", "complete"), "strain", "viral cRNA", False, 0, 8,
     "@@[organism]@@ @@[strain]@@ RNA, segment @@[segment]@@, complete sequence"),
    (_entry("segment", "4", "partial"), "strain", "viral cRNA", False, 0, 8,
     "@@[organism]@@ @@[strain]@@ RNA, segment @@[segment]@@"),
    # fallback (unknown type)
    (_entry("weird", "", "complete"), None, "genomic DNA", False, 0, 0,
     "@@[organism]@@ DNA, @@[entry]@@"),
    # mol token via mRNA
    (None, None, "mRNA", True, 0, 0,
     "@@[organism]@@ mRNA, @@[submitter_seqid]@@"),
    (None, None, "", True, 0, 0,
     "@@[organism]@@ DNA, @@[submitter_seqid]@@"),
])
def test_ff_definition_meta(entry, source_identifier, mol_type, is_wgs, chrom, seg, expected):
    out = ff_definition(entry, source_identifier, mol_type, is_wgs=is_wgs,
                        chromosome_count=chrom, segment_count=seg)
    assert out == expected


@pytest.mark.parametrize("entry,chrom,seg", [
    (_entry("plasmid", "", "complete"), 0, 0),
    (_entry("plasmid", "", "partial"), 0, 0),
    (_entry("chromosome", "", "complete"), 2, 0),
    (_entry("chromosome", "", "partial"), 2, 0),
    (_entry("segment", "", "complete"), 0, 2),
    (_entry("segment", "", "partial"), 0, 2),
])
def test_ff_definition_empty_seqname_raises(entry, chrom, seg):
    with pytest.raises(ValueError):
        ff_definition(entry, "strain", "genomic DNA", is_wgs=False,
                      chromosome_count=chrom, segment_count=seg)


def test_ff_definition_single_empty_seqname_allowed():
    # single chromosome / single segment with empty seq_name must NOT raise
    assert ff_definition(_entry("chromosome", "", "complete"), None, "genomic DNA",
                         chromosome_count=1) == "@@[organism]@@ DNA, chromosome, complete genome"
    assert ff_definition(_entry("segment", "", "complete"), None, "genomic DNA",
                         segment_count=1) == "@@[organism]@@ DNA, complete genome"
