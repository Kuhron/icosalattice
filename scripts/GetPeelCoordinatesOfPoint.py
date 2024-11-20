import math
import numpy as np
import random

import icosalattice.IcosahedronMath as icm
import icosalattice.MapCoordinateMath as mcm
import icosalattice.StartingPoints as sp
from icosalattice.UnitSpherePoint import UnitSpherePoint
import icosalattice.Faces as fc
import icosalattice.Edges as ed
from icosalattice.PointRepresentationAsFloat import point_code_to_float, point_float_to_code
import icosalattice.PeelCoordinates as pe


edge_midpoints = ed.get_edge_midpoints()


while True:
    # pick a test point
    if random.random() < 1/3:
        pc_orig = icm.get_random_point_code(min_iterations=0, expected_iterations=3, max_iterations=9)
        p = UnitSpherePoint.from_xyz(*pe.get_xyz_from_point_code_using_peel_coordinates(pc_orig))
        print(f"random point: {pc_orig}")
    elif random.random() < 1/2:
        i = random.randrange(len(sp.STARTING_POINT_CODES))
        pc_orig = sp.STARTING_POINT_CODES[i]
        p = sp.STARTING_POINTS[i]
        print(f"starting point {pc_orig}")
    else:
        edge_name, p = random.choice(list(edge_midpoints.items()))
        spc_expected, direction_expected = ed.get_ancestor_starting_point_and_direction_of_edge(edge_name)
        pc_orig = spc_expected + direction_expected
        print(f"point at middle of edge {edge_name}: {pc_orig}")
    
    spc, l, d = pe.get_peel_coordinates_of_point(p)
    print(f"{spc = }, {l = }, {d = } from p")
    pc = pe.get_point_code_from_peel_coordinates(spc, l, d)
    print(f"{pc = } from peel coords")
    assert pc == pc_orig

    ll = p.latlondeg(as_array=False)
    fs = fc.get_faces_of_point_by_closest_center(p)
    fpc = point_code_to_float(pc)
    print(f"{fpc = }")

    print(f"latlon {ll}\nis at peel coords L={l}, D={d} from point {spc}\ngot point code {pc}, float {fpc}")
    input("check")
    print()
