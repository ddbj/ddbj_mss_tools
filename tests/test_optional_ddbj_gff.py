"""ddbj-gff is an optional dependency needed only by gff2mss. Importing any CLI
(including gff2mss's) must not eagerly load ddbj_gff, so the other tools install
and run without it. Each check runs in a fresh subprocess to avoid import-cache
cross-talk, and is valid whether or not ddbj-gff happens to be installed.
"""

import os
import subprocess
import sys
from pathlib import Path

SRC = str(Path(__file__).resolve().parents[1] / "src")


def _import_check(imports: str):
    code = f"import sys; import {imports}; assert 'ddbj_gff' not in sys.modules; print('ok')"
    env = dict(os.environ)
    env["PYTHONPATH"] = SRC + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)


def test_non_gff_clis_import_without_ddbj_gff():
    r = _import_check("mss_builder.cli, egapx2mss.cli, mss2ff.cli, batch_wgs_builder.cli")
    assert r.returncode == 0, f"stdout={r.stdout!r} stderr={r.stderr!r}"


def test_gff2mss_package_import_does_not_eagerly_load_ddbj_gff():
    r = _import_check("gff2mss")
    assert r.returncode == 0, f"stdout={r.stdout!r} stderr={r.stderr!r}"


def test_gff2mss_cli_import_does_not_eagerly_load_ddbj_gff():
    r = _import_check("gff2mss.cli")
    assert r.returncode == 0, f"stdout={r.stdout!r} stderr={r.stderr!r}"
