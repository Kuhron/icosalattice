import pytest

import icosalattice.PointRepresentationAsFloat as pf


def test_point_code_to_float():
    test_cases = {
        "A": 0.0,
        "B": 1.0,
        "C": 2.0,
        "D": 3.0,
        "H2": 7 + 3/4,
        "F1": 5 + 2/4,
        "J3": 9 + 1/4,
        "K0001": 10 + 2/(4**4),
        "E231033": 4 + 3/4 + 1/(4**2) + 2/(4**3) + 0 + 1/(4**5) + 1/(4**6),
    }
    for pc_expected, fpc_expected in test_cases.items():
        fpc_got = pf.point_code_to_float(pc_expected)
        assert fpc_got == fpc_expected, f"expected fpc {fpc_expected} for pc {pc_expected} but got {fpc_got}"
        pc_got = pf.point_float_to_code(fpc_expected)
        assert pc_got == pc_expected, f"expected pc {pc_expected} for fpc {fpc_expected} but got {pc_got}"
