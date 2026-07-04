import os
from Bio import SeqIO
from ddbj_gff import parse
from ddbj_gff.normalize.normalize import normalize
from ddbj_gff.validate import validate
from gff2mss.convert import build_entry_features
from gff2mss.config import load_config

FIX = os.path.join(os.path.dirname(__file__), "mss_fixtures")


def _cfg(tmp_path):
    p = tmp_path / "cp.toml"
    p.write_text('[source]\norganism="Aliinostoc maniaoense"\nmol_type="genomic DNA"\n'
                 '[locus_tag]\nprefix="ACPZ3T"\n[cds]\ntransl_table=11\n'
                 '[transcript]\nemit_mrna=false\n', encoding="utf-8")
    cfg, _ = load_config(str(p))
    return cfg


def test_cp187952_origin_spanning_end_to_end(tmp_path):
    with open(os.path.join(FIX, "cp187952_origin.gff3")) as fh:
        doc = parse(fh.read())
    work, _ = normalize(doc)

    # canonicalization: origin-spanning modA no longer flagged out-of-region
    diags = validate(work)
    assert not any(d.code == "feature-outside-region" for d in diags)

    seqs = {rec.id: rec.seq for rec in
            SeqIO.parse(os.path.join(FIX, "cp187952_origin.fasta"), "fasta")}
    per_entry = build_entry_features(work, seqs, _cfg(tmp_path), [])
    feats = per_entry["CP187952.1"]

    # modA is emitted as a CDS (clean translation) at the wrapped location — not a misc_feature
    cds_locs = [f.location for f in feats if f.key == "CDS"]
    assert "complement(join(4447..5125,1..143))" in cds_locs
    assert not any(f.key == "misc_feature" for f in feats)
