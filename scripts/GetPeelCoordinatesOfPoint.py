import math
import numpy as np
import random

import icosalattice.IcosahedronMath as icm
import icosalattice.MapCoordinateMath as mcm
import icosalattice.StartingPoints as sp
from icosalattice.UnitSpherePoint import UnitSpherePoint
import icosalattice.Faces as fc
import icosalattice.Edges as ed


faces = fc.get_face_names()
starting_points = sp.STARTING_POINTS
edge_midpoints = ed.get_edge_midpoints()
labels = sp.STARTING_POINT_CODES
label_to_latlon = {label: p.latlondeg() for label, p in zip(labels, starting_points)}



def get_peel_coordinates_from_point_code(pc):
    raise NotImplementedError


def get_xyz_from_peel_coordinates(starting_pc, l_coord, d_coord):
    raise NotImplementedError


def get_peel_coordinates_from_xyz(xyz):
    # TODO move relevant stuff from get_peel_coordinates_of_point() to here
    raise NotImplementedError


def get_xyz_from_point_code_using_peel_coordinates(pc):
    spc, l, d = get_peel_coordinates_from_point_code(pc)
    xyz = get_xyz_from_peel_coordinates(spc, l, d)
    return xyz


def get_point_code_from_xyz_using_peel_coordinates(xyz):
    spc, l, d = get_peel_coordinates_from_xyz(xyz)
    pc = icm.get_point_code_from_peel_coordinates(spc, l, d)
    return pc



# optimizations I could make (but don't do it prematurely!)
# - for the vector composition stuff for getting peel coordinates from xyz and vice versa
# - - could be good to keep all the relevant displacement vectors as constants in one of the modules (maybe Edges)
# - - and ideally have their values analytically determined and specified like how MID_LAT_DEG is known from trig
# - - then can just look them up without recalculating them
# - rewrite the matrix equation code as more direct expressions for a1/a2/a3 rather than inverting a matrix (just write the inverse yourself)
# - - so that the math is done more directly with basic arithmetic operations in raw Python rather than NumPy


while True:
    # pick a test point
    if random.random() < 1/3:
        p = UnitSpherePoint.random()
        pc = None
        print(f"random point: {p}")
    elif random.random() < 1/2:
        i = random.randrange(len(sp.STARTING_POINT_CODES))
        pc = sp.STARTING_POINT_CODES[i]
        p = starting_points[i]
        print(f"point {pc}")
    else:
        edge_name, p = random.choice(list(edge_midpoints.items()))
        print(f"point at middle of edge {edge_name}")
    
    spc, l, d = icm.get_peel_coordinates_of_point(p)
    pc = icm.get_point_code_from_peel_coordinates(spc, l, d)
    ll = p.latlondeg(as_array=False)
    fs = fc.get_faces_of_point_by_closest_center(p)
    pc_bin = point_code_to_binary(pc)
    assert binary_to_point_code(pc_bin) == pc
    print(f"latlon {ll}\nis at peel coords L={l}, D={d} from point {spc}\ngot point code {pc}, binary {pc_bin}")
    input("check")
    print()
