import pytest

from egapx2mss import asn_tools

ASN = "Seq-entry ::= set {\n}\n"


def _fake_tool(tmp_path, name, script):
    p = tmp_path / name
    p.write_text("#!/bin/sh\n" + script)
    p.chmod(0o755)
    return p


@pytest.fixture
def no_download(monkeypatch):
    def fail(*a, **kw):
        raise AssertionError("must not re-download")
    monkeypatch.setattr(asn_tools, "_download_tool", fail)


@pytest.mark.parametrize("runner", ["gb", "fsa"])
def test_unrunnable_binary_reports_stderr(tmp_path, no_download, runner):
    # what the loader prints when a shared library is missing
    tool = _fake_tool(tmp_path, "tool", "echo 'error while loading shared libraries: libidn.so.11' >&2; exit 127\n")
    asn = tmp_path / "a.asn"
    asn.write_text(ASN)
    with pytest.raises(RuntimeError, match="libidn.so.11"):
        if runner == "gb":
            asn_tools.run_asn2gb_tbl(tool, str(asn), str(tmp_path / "a.tbl"))
        else:
            asn_tools.run_asn2fsa(tool, str(asn), str(tmp_path / "a.fa"), str(tmp_path))


def test_empty_output_error_includes_stderr(tmp_path, monkeypatch):
    tool = _fake_tool(tmp_path, "asn2gb", "echo 'This copy is too old' >&2; exit 1\n")
    monkeypatch.setattr(asn_tools, "_download_tool", lambda *a, **kw: tool)
    asn = tmp_path / "a.asn"
    asn.write_text(ASN)
    with pytest.raises(RuntimeError, match="too old"):
        asn_tools.run_asn2gb_tbl(tool, str(asn), str(tmp_path / "a.tbl"))
