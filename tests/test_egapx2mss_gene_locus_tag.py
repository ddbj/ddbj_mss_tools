from egapx2mss.ann_writer import write_ddbj_ann

TBL = (
    ">Feature lcl|chr1\n"
    "1\t30\tgene\n"
    "\t\t\tgene\tABC1\n"
    "\t\t\tgene_desc\tABC family member 1\n"
    "\t\t\tlocus_tag\tegapxtmp_000010\n"
    "1\t30\tmRNA\n"
    "\t\t\tgene\tABC1\n"
    "\t\t\tlocus_tag\tegapxtmp_000010\n"
    "\t\t\tproduct\tABC protein\n"
    "1\t30\tCDS\n"
    "\t\t\tgene\tABC1\n"
    "\t\t\tlocus_tag\tegapxtmp_000010\n"
    "\t\t\tproduct\tABC protein\n"
)
FASTA = ">chr1\n" + "ATG" * 10 + "\n"


def _ann(tmp_path, **kw):
    (tmp_path / "a.tbl").write_text(TBL)
    (tmp_path / "a.fa").write_text(FASTA)
    write_ddbj_ann(str(tmp_path / "a.tbl"), str(tmp_path / "a.fa"), str(tmp_path / "a.ann"), **kw)
    return (tmp_path / "a.ann").read_text()


def test_gene_omitted_by_default(tmp_path):
    ann = _ann(tmp_path)
    assert "\tgene\t1..30\t" not in ann
    assert "\tmRNA\t1..30\tlocus_tag\tegapxtmp_000010\n" in ann
    assert "\tCDS\t1..30\tlocus_tag\tegapxtmp_000010\n" in ann


def test_add_gene(tmp_path):
    assert "\tgene\t1..30\tlocus_tag\tegapxtmp_000010\n" in _ann(tmp_path, add_gene=True)


def test_locus_tag_prefix(tmp_path):
    ann = _ann(tmp_path, add_gene=True, locus_tag_prefix="NIGBMD")
    assert "egapxtmp" not in ann
    assert ann.count("locus_tag\tNIGBMD_000010\n") == 3


def test_gene_desc_moves_to_mrna_note(tmp_path):
    ann = _ann(tmp_path)
    assert ann.count("\tnote\tABC family member 1\n") == 1
    mrna = ann.split("\tmRNA\t")[1].split("\tCDS\t")[0]
    assert "\tnote\tABC family member 1\n" in mrna
    # with --add_gene the description stays on the gene only
    ann = _ann(tmp_path, add_gene=True)
    assert ann.count("\tnote\tABC family member 1\n") == 1
    assert "\tnote\t" not in ann.split("\tmRNA\t")[1]
