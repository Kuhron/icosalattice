import pytest
import numpy as np

import icosalattice.StartingPoints as sp
import icosalattice.Edges as ed
import icosalattice.PeelCoordinates as pe

from TestUtil import TEST_POINT_CODES


def test_peel_coordinates_vs_xyz():
    for pc0 in TEST_POINT_CODES:
        spc1_from_0, l1_from_0, d1_from_0 = pe.get_peel_coordinates_from_point_code(pc0)
        xyz1_from_1 = pe.get_xyz_from_peel_coordinates(spc1_from_0, l1_from_0, d1_from_0)
        spc1_from_1, l1_from_1, d1_from_1 = pe.get_peel_coordinates_from_xyz(xyz1_from_1)
        
        assert spc1_from_0 == spc1_from_1

        # potential float errors because of going through xyz
        assert np.isclose(l1_from_0, l1_from_1, atol=1e-9)
        assert np.isclose(d1_from_0, d1_from_1, atol=1e-9)

        pc1 = pe.get_point_code_from_peel_coordinates(spc1_from_0, l1_from_0, d1_from_0)
        assert pc1 == pc0
