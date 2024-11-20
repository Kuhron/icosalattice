import pytest
import numpy as np

import icosalattice.PeelCoordinates as pe
import icosalattice.IcosahedronMath as icm
import icosalattice.CoordinatesByAncestry as anc

from TestUtil import TEST_POINT_CODES


# TODO the math of position on the face plane vs on the actual sphere section for that face is wrong
# - need to account for the distortion introduced by warping the sphere section onto the face and vice versa
# - e.g. it is true that C1 is halfway between C and A on both the sphere section (by construction) and the face plane (because of symmetry of the sphere section)
# - and it is true that C11 is halfway between C1 and A on the sphere section (by construction)
# - but it is NOT true that C11 is halfway between C1 and A on the face plane! I think it is a little bit closer to C1


def test_xyz_by_ancestry_vs_peel_coordinates():
    for pc in TEST_POINT_CODES:
        print(f"testing {pc = }")
        xyz_by_ancestry = anc.get_xyz_from_point_code_using_ancestry(pc)
        xyz_by_peel_coordinates = pe.get_xyz_from_point_code_using_peel_coordinates(pc)
        print(f"{xyz_by_ancestry = }")
        print(f"{xyz_by_peel_coordinates = }")
        assert np.isclose(xyz_by_ancestry, xyz_by_peel_coordinates, atol=1e-9).all(), f"by ancestry: {xyz_by_ancestry}\nby peel coords: {xyz_by_peel_coordinates}"
        print()

