from gaopaw import parse_num_objs, parse_elems


def test_parse_num_objs():
    assert parse_num_objs('/scr/szymansk/gaopaw/tests/sample_workdir') == 15

def test_parse_elems():
    assert parse_elems('H2O2Na2FeCo34Ag5OSH2') == ['H','O','Na','Fe','Co','Ag','S']
