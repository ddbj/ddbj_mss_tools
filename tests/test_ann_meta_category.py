from ann_meta.category import CategoryGuess, EntrySketch, detect_category


def _entry(seq_id="seq1", has_topology=False, ff_definition=None, **qualifiers):
    return EntrySketch(seq_id=seq_id, has_topology=has_topology, ff_definition=ff_definition,
                       qualifiers=dict(qualifiers))


def _plain_entries():
    return [_entry(organism="Escherichia coli")]


def test_tls_by_datatype():
    g = detect_category({"DATATYPE": [("type", "TLS")]}, _plain_entries(), 500)
    assert g.value == "TLS"
    assert "DATATYPE type=TLS" in g.evidence
    assert g.tpa is False


def test_tls_keyword_beats_tsa_keyword():
    common = {"KEYWORD": [("keyword", "TSA; Targeted Locus Study")]}
    g = detect_category(common, _plain_entries(), 500)
    assert g.value == "TLS"


def test_tsa_by_division():
    g = detect_category({"DIVISION": [("division", "TSA")]}, _plain_entries(), 500)
    assert g.value == "TSA"
    assert "DIVISION division=TSA" in g.evidence


def test_tsa_by_keyword_case_insensitive():
    g = detect_category({"KEYWORD": [("keyword", "transcriptome shotgun assembly")]}, _plain_entries(), 500)
    assert g.value == "TSA"


def test_htg_by_keyword_phase():
    g = detect_category({"KEYWORD": [("keyword", "HTG"), ("keyword", "HTGS_PHASE1")]}, _plain_entries(), 500)
    assert g.value == "HTG"
    assert "KEYWORD HTGS_PHASE1" in g.evidence


def test_htg_by_keyword_phase_only():
    g = detect_category({"KEYWORD": [("keyword", "HTGS_PHASE2")]}, _plain_entries(), 500)
    assert g.value == "HTG"


def test_htc_by_keyword_fli():
    g = detect_category({"KEYWORD": [("keyword", "HTC"), ("keyword", "FLI_CDNA")]}, _plain_entries(), 500)
    assert g.value == "HTC"
    assert "KEYWORD FLI_CDNA" in g.evidence


def test_est_by_division():
    g = detect_category({"DIVISION": [("division", "EST")]}, _plain_entries(), 500)
    assert g.value == "EST"
    assert g.evidence == ["DIVISION division=EST"]


def test_mag_by_keyword_beats_wgs_datatype():
    common = {"DATATYPE": [("type", "WGS")],
              "KEYWORD": [("keyword", "WGS"), ("keyword", "STANDARD_DRAFT"), ("keyword", "MAG")]}
    g = detect_category(common, _plain_entries(), 500)
    assert g.value == "MAG"
    assert "KEYWORD MAG" in g.evidence


def test_mag_by_keyword_semicolon_split_case_insensitive():
    common = {"KEYWORD": [("keyword", "STANDARD_DRAFT; metagenome assembled genome")]}
    g = detect_category(common, _plain_entries(), 500)
    assert g.value == "MAG"


def test_mag_by_env_division_and_metagenome_source():
    common = {"DATATYPE": [("type", "WGS")], "DIVISION": [("division", "ENV")]}
    entries = [_entry(organism="uncultured Bacteroides sp.", metagenome_source="human gut metagenome")]
    g = detect_category(common, entries, 500)
    assert g.value == "MAG"
    assert g.evidence == ["DIVISION division=ENV", "source metagenome_source"]


def test_env_division_without_metagenome_source_is_not_mag():
    common = {"DATATYPE": [("type", "WGS")], "DIVISION": [("division", "ENV")]}
    g = detect_category(common, _plain_entries(), 500)
    assert g.value == "WGS"


def test_sag_by_keyword():
    common = {"DATATYPE": [("type", "WGS")], "KEYWORD": [("keyword", "Single Amplified Genome")]}
    g = detect_category(common, _plain_entries(), 500)
    assert g.value == "SAG"
    assert "KEYWORD Single Amplified Genome" in g.evidence


def test_wgs_by_datatype():
    g = detect_category({"DATATYPE": [("type", "WGS")]}, _plain_entries(), 500)
    assert g.value == "WGS"
    assert g.evidence == ["DATATYPE type=WGS"]


def test_wgs_by_keyword_only():
    g = detect_category({"KEYWORD": [("keyword", "wgs")]}, _plain_entries(), 500)
    assert g.value == "WGS"
    assert g.evidence == ["KEYWORD wgs"]


def test_gnm_by_topology():
    entries = [_entry(has_topology=True, organism="Escherichia coli")]
    g = detect_category({}, entries, 500)
    assert g.value == "GNM"
    assert g.evidence == ["TOPOLOGY present"]


def test_gnm_by_length():
    g = detect_category({}, _plain_entries(), 1_000_001)
    assert g.value == "GNM"
    assert g.evidence == ["first entry 1,000,001 bp > 1,000,000"]


def test_gnm_not_by_length_at_threshold():
    g = detect_category({}, _plain_entries(), 1_000_000)
    assert g.value is None


def test_gnm_by_ff_definition_complete_genome():
    entries = [_entry(ff_definition="@@[organism]@@ @@[strain]@@ DNA, complete genome", organism="E. coli")]
    g = detect_category({}, entries, 500)
    assert g.value == "GNM"
    assert g.evidence == ["ff_definition 'complete genome'"]


def test_gnm_by_ff_definition_chromosome():
    entries = [_entry(ff_definition="Escherichia coli DNA, Chromosome 1", organism="E. coli")]
    g = detect_category({}, entries, 500)
    assert g.value == "GNM"
    assert g.evidence == ["ff_definition ', chromosome'"]


def test_gnm_chromosome_phrase_is_anchored():
    # a gene named after chromosomes is not a chromosome entry
    entries = [_entry(ff_definition="@@[organism]@@ chromosome partition protein gene, partial cds",
                     organism="E. coli")]
    assert detect_category({}, entries, 500).value is None


def test_gnm_requires_no_datatype():
    entries = [_entry(has_topology=True, organism="Escherichia coli")]
    g = detect_category({"DATATYPE": [("type", "WGS")]}, entries, 5_000_000)
    assert g.value == "WGS"


def test_misc_plasmid_only():
    entries = [
        _entry("p1", ff_definition="Escherichia coli plasmid pA DNA, complete sequence", organism="E. coli"),
        _entry("p2", ff_definition="Escherichia coli plasmid pB DNA, complete sequence", organism="E. coli"),
    ]
    g = detect_category({}, entries, 50_000)
    assert g.value == "MISC"
    assert g.evidence == ["ff_definition 'plasmid ... complete sequence' (all entries)"]


def test_misc_plasmid_not_all_entries():
    entries = [
        _entry("p1", ff_definition="Escherichia coli plasmid pA DNA, complete sequence", organism="E. coli"),
        _entry("p2", ff_definition="Escherichia coli gene X, partial cds", organism="E. coli"),
    ]
    g = detect_category({}, entries, 50_000)
    assert g.value is None


def test_misc_organelle_only():
    entries = [_entry("m1", organism="Homo sapiens", organelle="mitochondrion")]
    g = detect_category({}, entries, 16_569)
    assert g.value == "MISC"
    assert g.evidence == ["source organelle (all entries)"]


def test_misc_virus_organism():
    entries = [_entry(organism="Tobacco mosaic virus")]
    g = detect_category({}, entries, 6_395)
    assert g.value == "MISC"
    assert g.evidence == ["organism contains 'virus'"]


def test_misc_phage_organism():
    entries = [_entry(organism="Escherichia Phage T4")]
    g = detect_category({}, entries, 5_000)
    assert g.value == "MISC"
    assert g.evidence == ["organism contains 'phage'"]


def test_none_for_marker_set():
    entries = [_entry(ff_definition="Lepidoptera sp. COI gene for cytochrome c oxidase subunit I, partial cds",
                      organism="Lepidoptera sp.")]
    g = detect_category({}, entries, 500)
    assert g.value is None
    assert g.evidence == []
    assert g.tpa is False


def test_tpa_flag():
    common = {"DATATYPE": [("type", "TPA-WGS")], "KEYWORD": [("keyword", "TPA"), ("keyword", "MAG")]}
    g = detect_category(common, _plain_entries(), 500)
    assert isinstance(g, CategoryGuess)
    assert g.tpa is True
    assert g.value is None
    assert g.evidence == ["DATATYPE type=TPA-WGS"]


def test_misc_organelle_only_beats_gnm_topology():
    entries = [_entry("m1", has_topology=True, ff_definition="Homo sapiens mitochondrial DNA, complete genome",
                      organism="Homo sapiens", organelle="mitochondrion")]
    g = detect_category({}, entries, 16_569)
    assert g.value == "MISC"
    assert g.evidence == ["source organelle (all entries)"]


def test_misc_virus_beats_gnm_complete_genome():
    entries = [_entry(ff_definition="Tomato mosaic virus RNA, complete genome", organism="Tomato mosaic virus")]
    g = detect_category({}, entries, 6_400)
    assert g.value == "MISC"
    assert g.evidence == ["organism contains 'virus'"]


def test_misc_plasmid_only_with_topology():
    entries = [_entry("p1", has_topology=True, ff_definition="Escherichia coli plasmid pX DNA, complete sequence",
                      organism="Escherichia coli")]
    g = detect_category({}, entries, 60_000)
    assert g.value == "MISC"
    assert g.evidence == ["ff_definition 'plasmid ... complete sequence' (all entries)"]
