import json
import sys
from pathlib import Path

import pytest

from tsa2mss.cli import main

EX = Path(__file__).parent.parent / "examples" / "tsa2mss"


def run(monkeypatch, tmp_path, common, *extra):
    monkeypatch.setattr(sys, "argv", ["tsa2mss", str(EX / "transcripts.fa"), "--common", str(common),
                                      "-o", str(tmp_path), "-p", "out", *extra])
    main()
    return (tmp_path / "out.ann").read_text()


def test_tsa(monkeypatch, tmp_path):
    ann = run(monkeypatch, tmp_path, EX / "common_tsa.json")
    assert "\tDIVISION\t\tdivision\tTSA\n" in ann
    assert "\tKEYWORD\t\tkeyword\tTSA\n\t\t\tkeyword\tTranscriptome Shotgun Assembly\n" in ann
    assert "\tsource\t1..E\tmol_type\ttranscribed RNA\n" in ann
    assert "ff_definition\t@@[organism]@@ @@[strain]@@ RNA, @@[submitter_seqid]@@\n" in ann
    assert "ENV" not in ann and "environmental_sample" not in ann
    assert "contig2\tassembly_gap\t22..31\testimated_length\tknown\n" in ann
    assert (tmp_path / "out.fa").read_text().count("//") == 2


def test_meta_tsa(monkeypatch, tmp_path):
    ann = run(monkeypatch, tmp_path, EX / "common_meta_tsa.json", "--meta")
    assert "\tDIVISION\t\tdivision\tTSA\n" in ann
    assert "\t\t\tkeyword\tENV\n" in ann
    assert "\t\t\tenvironmental_sample\t\n" in ann
    assert "\t\t\tisolation_source\tcoral reef\n" in ann


def test_missing_required_filled_empty(monkeypatch, tmp_path, capsys):
    common = json.loads((EX / "common_meta_tsa.json").read_text())
    del common["SOURCE"]["isolation_source"], common["SOURCE"]["collection_date"]
    del common["DBLINK"]["sequence read archive"]
    p = tmp_path / "c.json"
    p.write_text(json.dumps(common))
    ann = run(monkeypatch, tmp_path, p, "--meta")
    err = capsys.readouterr().err
    for key in ("isolation_source", "collection_date", "sequence read archive"):
        assert f"'{key}' is missing" in err
        assert f"\t{key}\t\n" in ann


def test_unknown_gap_length_rejected(monkeypatch, tmp_path):
    common = json.loads((EX / "common_tsa.json").read_text())
    common["ASSEMBLY_GAP"]["linkage_evidence"] = "align genus"
    p = tmp_path / "c.json"
    p.write_text(json.dumps(common))
    with pytest.raises(SystemExit, match="estimated_length 'unknown'"):
        run(monkeypatch, tmp_path, p)


def test_st_comment_tagset(monkeypatch, tmp_path, capsys):
    common = json.loads((EX / "common_tsa.json").read_text())
    st = common.pop("ST_COMMENT")
    p = tmp_path / "c.json"
    p.write_text(json.dumps(common))
    ann = run(monkeypatch, tmp_path, p)
    assert "\tST_COMMENT\t\ttagset_id\tAssembly-Data\n" in ann

    common["ST_COMMENT"] = {**st, "tagset_id": "Genome-Assembly-Data"}
    p.write_text(json.dumps(common))
    capsys.readouterr()
    run(monkeypatch, tmp_path, p)
    assert "tagset_id should be 'Assembly-Data'" in capsys.readouterr().err
