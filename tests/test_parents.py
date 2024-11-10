import pytest

from icosalattice.ParentsAndChildren import get_parent_from_point_code


def test_parents():
    point_to_par = {
        "A": None,
        "J": None,
        "C0": "C",
        "C1": "C",
        "C1323": "C132",
        "C22": "C2",
        "C32": "C3",
        "C323": "C32",
        "D132": "D13",
        "D2132": "D213",
        "F221": "F22",
        "F3221": "F322",
        "G12": "G1",
        "G122": "G12",
        "G1221": "G122",
        "G2": "G",
        "G22": "G2",
        "G221": "G22",
        "K23020": "K2302",
    }
    for pc, par in point_to_par.items():
        par2 = get_parent_from_point_code(pc)
        assert par == par2, f"par of {pc} is {par} but got {par2}"
