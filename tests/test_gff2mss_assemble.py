from gff2mss.assemble import build_ann_text


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return str(p)


def test_organelle_circular_source_and_topology(tmp_path):
    gff = _write(tmp_path, "o.gff3",
        "##gff-version 3\n"
        "CP\tLiftoff\tgene\t1\t9\t.\t+\t.\tID=g1;gene=rpl5\n"
        "CP\tLiftoff\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        "CP\tLiftoff\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g1.t1;product=50S ribosomal protein L5;transl_table=11\n")
    fasta = _write(tmp_path, "o.fa", ">CP\nATGAAATAA\n")
    mss_cfg = _write(tmp_path, "o.toml",
        '[source]\norganism="x"\nmol_type="genomic DNA"\n[locus_tag]\nprefix="HAKA"\n[cds]\ntransl_table=1\n')
    common = _write(tmp_path, "c.json",
        '{"DATE":{"hold_date":"2099-12-31"},'
        '"SOURCE":{"organism":"Heterosigma akashiwo","mol_type":"genomic DNA"},'
        '"SOURCE_IDENTIFIER":"strain"}')
    roles = _write(tmp_path, "r.tsv", "CP\torganelle\tplastid:chloroplast\tcomplete\tcircular\n")
    text, _ = build_ann_text(gff, fasta, mss_cfg, common, roles, "GNM")
    assert "COMMON" in text
    assert "TOPOLOGY" in text and "circular" in text
    assert "\torganelle\tplastid:chloroplast" in text
    assert "50S ribosomal protein L5" in text
    assert "\tCDS\t" in text


def test_nuclear_wgs_submitter_seqid(tmp_path):
    gff = _write(tmp_path, "n.gff3",
        "##gff-version 3\n"
        "scaffold_1\tS\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "scaffold_1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        "scaffold_1\tS\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g1.t1\n")
    fasta = _write(tmp_path, "n.fa", ">scaffold_1\nATGAAATAA\n")
    mss_cfg = _write(tmp_path, "n.toml",
        '[source]\norganism="x"\nmol_type="genomic DNA"\n[locus_tag]\nprefix="HAKA"\n[cds]\ntransl_table=1\n[product]\ndefault="hypothetical protein"\n')
    common = _write(tmp_path, "c.json",
        '{"DATE":{"hold_date":"2099-12-31"},'
        '"SOURCE":{"organism":"Heterosigma akashiwo","mol_type":"genomic DNA"},'
        '"SOURCE_IDENTIFIER":"strain"}')
    # no roles file -> all unplaced -> WGS
    text, _ = build_ann_text(gff, fasta, mss_cfg, common, None, "WGS")
    assert "submitter_seqid\tscaffold_1" in text
    assert "TOPOLOGY" not in text          # linear WGS
    assert "hypothetical protein" in text  # no product -> default
