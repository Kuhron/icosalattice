import pytest
import numpy as np

import icosalattice.StartingPoints as sp
import icosalattice.Edges as ed
import icosalattice.PeelCoordinates as pe
from icosalattice.MapCoordinateMath import get_latlon_from_xyz, get_xyz_from_latlon

from TestUtil import TEST_POINT_CODES, get_test_floats_01


def test_invert_ld_adjustment():
    n_points = 100
    floats = get_test_floats_01(2 * n_points)

    for i in range(n_points):
        spc1 = "C"
        l1, d1, *floats = floats
        f12 = pe.get_adjusted_peel_coordinates_from_raw_peel_coordinates
        f21 = pe.get_raw_peel_coordinates_from_adjusted_peel_coordinates  # FIXME BUG HERE
        spc12, l12, d12 = f12((spc1, l1, d1))
        spc121, l121, d121 = f21((spc1, l12, d12))
        # assert spc1 == spc12 == spc121
        # assert np.isclose(l1, l121, rtol=1e-12)
        # assert np.isclose(d1, d121, rtol=1e-12)
        print(f"\nl: {l1:.6f} {l12:.6f} {l121:.6f}")
        print(f"d: {d1:.6f} {d12:.6f} {d121:.6f}")
        input("check\n")
    
    assert len(floats) == 0, "didn't use up all the floats generated"
