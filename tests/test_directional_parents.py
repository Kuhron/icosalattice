import pytest

from icosalattice.ParentsAndChildren import get_directional_parent_from_point_code


def test_directional_parents():
    point_to_dpar = {
        "C0": "C",
        "C1": "A",
        "C1323": "C201",
        "C22": "K",
        "C32": "L1",
        "C323": "L01",
        "D132": "D21",
        "D2132": "D221",
        "F221": "E33",
        "F3221": "F233",
        "G12": "E1",
        "G122": "E1",
        "G1221": "E101",
        "G2": "E",
        "G22": "E",
        "G221": "E01",
        "K23020": "K2302",
    }
    for pc, dpc in point_to_dpar.items():
        dpc2 = get_directional_parent_from_point_code(pc)
        assert dpc == dpc2, f"dpar of {pc} is {dpc} but got {dpc2}"
