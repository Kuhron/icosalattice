import pytest


def test_directional_parents():
    point_to_dpar = {
        "D2132": "D221",
        "D132": "D21",
        "G1221": "E101",
        "G221": "E01",
        "C1323": "C201",
        "C323": "L01",
        "F3221": "F233",
        "F221": "E33",
    }
    for pc, dpc in point_to_dpar.items():
        dpc2 = get_directional_parent(pc)
        assert dpc == dpc2
