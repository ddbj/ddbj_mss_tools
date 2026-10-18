import json
import shutil
from pathlib import Path

from ann_fix.annfile import AnnFile
from ann_fix.cli import main
from ann_fix.report import Finding

FX = Path(__file__).parent / "ann_fix_fixtures"


def _report_for(ann_path: Path) -> str:
    """A details report whose Line[N] values point at the mini.ann features."""
    ann = AnnFile.read(ann_path)
    cds2 = [f for f in ann.features() if f.key == "CDS"][1]
    trna = next(f for f in ann.features() if f.key == "tRNA")
    return "\n".join([
        "=== Validation Details (1 file set) ===",
        f"ANN5320:ERR:NC_001133.9:tRNA:Line[{trna.first + 1}]:The tRNA feature must not be completely contained within CDS exons. (Found: tRNA 600..671 in CDS 538..792)",
        f"ANN4300:WAR:NC_001133.9:CDS:Line[{cds2.first + 1}]:Complete CDS must be at least 90 bases (30 aa) long unless supported by experiment or inference. (Length: 60 bases)",
        f"AXS6810:WAR:NC_001133.9:CDS:Line[{cds2.first + 1}]:Non-canonical splice sites: GT-AG/GC-AG rule violation. (Found: GT-AT in CDS at 538..792)",
        "ANN2542:WAR:ALL:ALL:locus_tag:[FILE]:Duplicate locus_tag across features without common gene qualifiers. (locus_tag: 'SCEX_000010', genes: none)",
        "ANN0470:WAR:ALL:DBLINK:[FILE]:Missing SRA Run accession.",
    ]) + "\n"


def _setup(tmp_path):
    ann = tmp_path / "x.ann"
    shutil.copy(FX / "mini.ann", ann)
    rep = tmp_path / "details.txt"
    rep.write_text(_report_for(ann))
    return ann, rep


def test_detect_json(tmp_path, capsys):
    ann, rep = _setup(tmp_path)
    assert main(["--ann", str(ann), "--report", str(rep), "--detect", "--json"]) == 0
    data = {d["code"]: d for d in json.loads(capsys.readouterr().out)}
    assert data["ANN5320"]["default"] == "overlap-drop" and data["ANN5320"]["level"] == "ERR"
    assert data["ANN5320"]["count"] == 1 and data["ANN5320"]["examples"][0]["key"] == "tRNA"
    assert data["ANN5320"]["examples"][0]["locus_tag"] == "SCEX_000030"
    assert data["ANN2542"]["default"].startswith("ignore")
    assert data["ANN0470"]["default"] == "none"


def test_detect_text_lists_examples(tmp_path, capsys):
    ann, rep = _setup(tmp_path)
    assert main(["--ann", str(ann), "--report", str(rep)]) == 0
    out = capsys.readouterr().out
    assert "ANN4300" in out and "n=1" in out and "default=misc" in out and "538..792" in out


def test_apply_defaults_then_idempotent(tmp_path, capsys):
    ann, rep = _setup(tmp_path)
    out = tmp_path / "o.ann"
    assert main(["--ann", str(ann), "--report", str(rep), "--apply", "defaults", "--out", str(out)]) == 0
    printed = capsys.readouterr().out
    assert "re-validate:" in printed
    fixed = AnnFile.read(out)
    keys = [f.key for f in fixed.features() if f.entry != "COMMON"]
    assert keys == ["source", "mRNA", "CDS", "misc_feature"]       # tRNA dropped, 2nd pair -> misc
    misc = next(f for f in fixed.features() if f.key == "misc_feature")
    assert ("artificial_location", "low-quality sequence region") not in misc.qualifiers(fixed.rows)
    # applying the same report again on the fixed file: the report no longer matches -> nothing written, exit 1
    before = out.read_bytes()
    assert main(["--ann", str(out), "--report", str(rep), "--apply", "defaults", "--out", str(out)]) == 1
    assert out.read_bytes() == before
    assert "does not match" in capsys.readouterr().err


def test_apply_explicit_spec(tmp_path):
    ann, rep = _setup(tmp_path)
    out = tmp_path / "o.ann"
    assert main(["--ann", str(ann), "--report", str(rep),
                 "--apply", "ANN5320=overlap-drop:other,AXS6810=artificial-location:heterogeneous population sequenced",
                 "--out", str(out)]) == 0
    fixed = AnnFile.read(out)
    keys = [f.key for f in fixed.features() if f.entry != "COMMON"]
    assert keys == ["source", "mRNA", "CDS", "tRNA"]
    # AXS6810 pointed at the CDS that overlap-drop removed in the same run -> skipped, no artificial_location anywhere
    assert not any(k == "artificial_location" for f in fixed.features() for k, _ in f.qualifiers(fixed.rows))


def test_bad_keep_or_value_exits_2(tmp_path):
    ann, rep = _setup(tmp_path)
    out = str(tmp_path / "o.ann")
    assert main(["--ann", str(ann), "--report", str(rep), "--apply", "ANN5320=overlap-drop:bogus", "--out", out]) == 2
    assert main(["--ann", str(ann), "--report", str(rep), "--apply", "AXS6810=artificial-location:made up", "--out", out]) == 2
    assert main(["--ann", str(ann), "--report", str(rep), "--apply", "ANN5320=overlap-drop:CDS", "--out", out]) == 0


def test_apply_unknown_code_is_notice(tmp_path, capsys):
    ann, rep = _setup(tmp_path)
    out = tmp_path / "o.ann"
    assert main(["--ann", str(ann), "--report", str(rep), "--apply", "ANN9999=misc", "--out", str(out)]) == 0
    assert "no findings for ANN9999" in capsys.readouterr().out
    assert out.read_bytes() == ann.read_bytes()


def test_bad_operation_exits_2(tmp_path):
    ann, rep = _setup(tmp_path)
    assert main(["--ann", str(ann), "--report", str(rep), "--apply", "ANN4300=bogus", "--out", str(tmp_path / "o.ann")]) == 2
    assert main(["--ann", str(ann), "--report", str(rep), "--apply", "ANN4300=misc"]) == 2   # --out required
