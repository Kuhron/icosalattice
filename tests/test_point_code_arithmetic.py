import pytest

import icosalattice.PointCodeArithmetic as pca


def test_point_code_arithmetic():
    test_cases = {
        "C": {1: "A", 2: "K", 3: "L", -1: "D", -2: "E", -3: None},
        "D": {1: "C", 2: "L", 3: "B", -1: None, -2: "F", -3: "E"},
        "G": {1: "A", 2: "E", 3: "F", -1: "H", -2: "I", -3: None},
        "C1": {1: "A0", 2: "K1", 3: "C2", -1: "C0", -2: "E2", -3: "E1"},
        "D3": {1: "D2", 2: "L3", 3: "B0", -1: "F3", -2: "F2", -3: "D0"},
        "C100": {-2: "E211", -3: "E122"},
        "D333": {-1: "F333", -2: "F332"},
        "C110": {-2: "E121", -3: "E112"},
        "D033": {-1: "F233", -2: "F232"},
        "C011": {-2: "E212", -3: "E211"},
        "D300": {-1: "F322", -2: "F233"},
        "K121": {1: "I110", 2: "I101", 3: "K122"},
        "K2": {3: "J1"},
        "K22": {3: "J11"},
        "K23": {2: "J11", 3: "J10"},
    }
    for pc, d in test_cases.items():
        for x, target in d.items():
            got = pca.add_direction_to_point_code(pc, x)
            # print(f"{pc} {x:+} = {got}")
            if got != target:
                raise Exception(f"{pc} {x:+} = {target} but got {got}")
