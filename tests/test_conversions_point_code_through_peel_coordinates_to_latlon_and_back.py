import pytest
import numpy as np

import icosalattice.StartingPoints as sp
import icosalattice.Edges as ed
import icosalattice.PeelCoordinates as pe
from icosalattice.MapCoordinateMath import get_latlon_from_xyz, get_xyz_from_latlon

from TestUtil import TEST_POINT_CODES

# test that the inverses of conversions work correctly
# pc -> ld_raw -> pc
# ld_raw -> ld_adjusted -> ld_raw
# ld_adjusted -> xyz -> ld_adjusted
# xyz -> latlon -> xyz
# and finally, compose them all: pc -> ld_raw -> ld_adjusted -> xyz -> latlon and then reverse them all and check near-equality along each return step

# use path notation with numbers for each of the steps it's been to:
# 1 = point code
# 2 = raw peel coordinates
# 3 = corrected peel coordinates
# 4 = xyz
# 5 = latlon

def test_conversions_point_code_through_peel_coordinates_to_latlon_and_back():
    f12 = pe.get_raw_peel_coordinates_from_point_code
    f21 = lambda sld: pe.get_point_code_from_raw_peel_coordinates(sld, max_iterations=8, allow_clipping=True)

    f23 = pe.get_adjusted_peel_coordinates_from_raw_peel_coordinates
    f32 = pe.get_raw_peel_coordinates_from_adjusted_peel_coordinates

    f34 = pe.get_xyz_from_adjusted_peel_coordinates
    f43 = pe.get_adjusted_peel_coordinates_from_xyz

    f45 = get_latlon_from_xyz
    f54 = get_xyz_from_latlon

    for pc_1 in TEST_POINT_CODES:
        print(f"converting for pc = {pc_1}")

        sld_raw_12 = f12(pc_1)
        print(f"{sld_raw_12 = }")
        pc_121 = f21(sld_raw_12)
        print(f"{pc_121 = }")

        sld_adjusted_123 = f23(sld_raw_12)
        print(f"{sld_adjusted_123 = }")
        sld_raw_1232 = f32(sld_adjusted_123)
        print(f"{sld_raw_1232 = }")
        pc_12321 = f21(sld_raw_1232)
        print(f"{pc_12321 = }")

        xyz_1234 = f34(sld_adjusted_123)
        print(f"{xyz_1234 = }")
        sld_adjusted_12343 = f43(xyz_1234)
        print(f"{sld_adjusted_12343 = }")
        sld_raw_123432 = f32(sld_adjusted_12343)
        print(f"{sld_raw_123432 = }")
        pc_1234321 = f21(sld_raw_123432)
        print(f"{pc_1234321 = }")

        latlon_12345 = f45(xyz_1234)
        print(f"{latlon_12345 = }")
        xyz_123454 = f54(latlon_12345)
        print(f"{xyz_123454 = }")
        sld_adjusted_1234543 = f43(xyz_123454)
        print(f"{sld_adjusted_1234543 = }")
        sld_raw_12345432 = f32(sld_adjusted_1234543)
        print(f"{sld_raw_12345432 = }")
        pc_123454321 = f21(sld_raw_12345432)
        print(f"{pc_123454321 = }")

        pcs_to_check = [pc_1, pc_121, pc_12321, pc_1234321, pc_123454321]
        slds_raw_to_check = [sld_raw_12, sld_raw_1232, sld_raw_123432, sld_raw_12345432]
        slds_adjusted_to_check = [sld_adjusted_123, sld_adjusted_12343, sld_adjusted_1234543]
        xyzs_to_check = [xyz_1234, xyz_123454]

        print("\ndebug checking equality after all conversions:")
        print("pcs_to_check:")
        for x in pcs_to_check:
            print(f"\t{x}")
        print("slds_raw_to_check:")
        for x in slds_raw_to_check:
            print(f"\t{x}")
        print("slds_adjusted_to_check:")
        for x in slds_adjusted_to_check:
            print(f"\t{x}")
        print("xyzs_to_check:")
        for x in xyzs_to_check:
            print(f"\t{x}")
        print("\n--------\n")

        pc_ref = pcs_to_check[0]
        for pc in pcs_to_check[1:]:
            assert pc == pc_ref, pc
        spc_ref, l_raw_ref, d_raw_ref = slds_raw_to_check[0]
        for spc, l, d in slds_raw_to_check[1:]:
            assert spc == spc_ref, spc
            assert np.isclose(l, l_raw_ref, atol=1e-9), l
            assert np.isclose(d, d_raw_ref, atol=1e-9), d
        spc_ref, l_adjusted_ref, d_adjusted_ref = slds_adjusted_to_check[0]
        for spc, l, d in slds_adjusted_to_check[1:]:
            assert spc == spc_ref, spc
            assert np.isclose(l, l_adjusted_ref, atol=1e-9), l
            assert np.isclose(d, d_adjusted_ref, atol=1e-9), d
        assert all(np.isclose(xyz, xyzs_to_check[0], atol=1e-9).all() for xyz in xyzs_to_check[1:]), xyzs_to_check
