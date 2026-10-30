from common.qualifiers import known_qualifier_names, qualifier_def


def test_known_names_has_organism_and_44_entries():
    names = known_qualifier_names()
    assert isinstance(names, frozenset)
    assert "organism" in names
    assert len(names) == 44


def test_qualifier_def_flags():
    assert qualifier_def("environmental_sample")["boolean"] is True
    assert qualifier_def("strain")["name"] == "strain"
    assert qualifier_def("nope") is None
