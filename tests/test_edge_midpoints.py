import math
import numpy as np

import icosalattice.IcosahedronMath as icm
import icosalattice.MapCoordinateMath as mcm
import icosalattice.StartingPoints as sp
from icosalattice.Edges import get_edge_midpoints


def test_edge_midpoints():
    starting_points, adj = sp.STARTING_POINTS_AND_ADJACENCY
    labels = sp.STARTING_POINT_CODES
    label_to_latlon = {label: p.latlondeg() for label, p in zip(labels, starting_points)}

    midpoints_latlon = get_edge_midpoints()

    for (pc0, pc1), p in midpoints_latlon.items():
        mp_latlon_expected = p.latlondeg()
        ll0 = label_to_latlon[pc0]
        ll1 = label_to_latlon[pc1]
        mp_latlon_got = mcm.get_unit_sphere_midpoint_from_latlon(ll0, ll1, as_array=False)
        assert mp_latlon_got == mcm.get_unit_sphere_midpoint_from_latlon(ll1, ll0, as_array=False)
        # print(pc0, ll0, f"{mp_latlon_expected} vs {mp_latlon_got}", ll1, pc1)
        assert abs(mp_latlon_expected[0] - mp_latlon_got[0]) < 1e-9
        assert abs(mp_latlon_expected[1] - mp_latlon_got[1]) < 1e-9
