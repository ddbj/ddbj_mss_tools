"""End-to-end: type=segment in sequence roles -> mss_builder .ann output."""

from common.models import CommonModel
from common.source_builder import SequenceRoleEntry
from mss_builder.ann_writer import write_mss_ann


def _common():
    return CommonModel.model_validate({
        "DBLINK": {"project": "P", "sample": "S"},
        "SOURCE": {"organism": "Influenza A virus", "mol_type": "viral cRNA"},
    })


def test_multi_segment_emits_segment_qualifier(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">seg4\nACGTACGTAC\n>seg6\nTTTTGGGGCC\n")
    ann = tmp_path / "out.ann"
    roles = {
        "seg4": SequenceRoleEntry("seg4", "segment", "4", "complete", False),
        "seg6": SequenceRoleEntry("seg6", "segment", "6", "complete", False),
    }
    write_mss_ann(str(fasta), str(ann), common=_common(), sequence_roles=roles)
    text = ann.read_text()
    assert "\t\t\tsegment\t4\n" in text
    assert "Influenza A virus RNA, segment 4, complete sequence" in text


def test_single_segment_omits_qualifier_and_uses_genome(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">seg1\nACGTACGTAC\n")
    ann = tmp_path / "out.ann"
    roles = {"seg1": SequenceRoleEntry("seg1", "segment", "", "complete", False)}
    write_mss_ann(str(fasta), str(ann), common=_common(), sequence_roles=roles)
    text = ann.read_text()
    assert "Influenza A virus RNA, complete genome" in text
    assert "\t\t\tsegment\t" not in text
