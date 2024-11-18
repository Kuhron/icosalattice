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


def get_peel_coordinates_of_point(p):
    # peel coordinates are within the shape formed by combining this face and its up-facing or down/facing counterpart
    # and shearing those such that it's a square with the 0 point (C,D,E,F,G,H,I,J,K,L) on the upper right
    # the 1 direction going left, the 2 direction going down-left, and the 3 direction going down
    # so each peel has a north and a south square, each the parental watershed of one of the 10 mid-latitude starting points
    # within this square, measure how far left we go and how far down we go (each in interval from 0 to 1)
    # from this we should be able to get the point code directly, I hope?

    p_xyz = p.xyz(as_array=True)

    fs = fc.get_faces_of_point_by_closest_center(p)
    assert 1 <= len(fs) <= 5, f"got {len(fs)} faces for point, which should be impossible; {p = }"

    if len(fs) == 1:
        face ,= fs
        spc = face[0]
        l_coord, d_coord = get_peel_coordinates_of_point_from_face_name(p_xyz, face)
    else:
        # A and B (poles) are not on any peel, should return some special value to indicate this
        # starting points are at (0, 0) on their own peel
        # edge points in 1 direction from a starting point are on that peel (e.g. C1 is at D=0 from C)
        # edge points in 2 direction from a starting point should have the same peel coordinates whether you choose the upward or downward face (e.g. C2)
        # edge points in 3 direction from a starting point are on that peel (e.g. D3 is at L=0 from D)
        # edge points on the bottom or left of a peel are NOT on that peel, they're part of the one west of it
        
        if len(fs) == 5:
            # starting point
            if all("A" in x for x in fs):
                spc = "A"
            elif all("B" in x for x in fs):
                spc = "B"
            else:
                s = get_vertices_in_common_to_faces(fs)
                assert len(s) == 1, s
                spc ,= s
            l_coord = 0
            d_coord = 0
        elif len(fs) == 2:
            s = get_vertices_in_common_to_faces(fs)
            spc, direction = ed.get_ancestor_starting_point_and_direction_of_edge(s)
            # print(f"{spc = }, {direction = }")

            if direction == "1":
                # get peel coordinates on up-pointing face radiating from ancestor point
                face ,= [x for x in fs if x[0] == spc]
                l_coord, d_coord = get_peel_coordinates_of_point_from_face_name(p_xyz, face)
            elif direction == "3":
                # edge point in the 2 direction, so it borders both the up-pointing and down-pointing faces
                # verify that getting coords from both faces give same result
                face ,= [x for x in fs if x[0] == spc]
                l_coord, d_coord = get_peel_coordinates_of_point_from_face_name(p_xyz, face)
            elif direction == "2":
                # get peel coordinates on down-pointing face radiating from ancestor point
                l_coord0, d_coord0 = get_peel_coordinates_of_point_from_face_name(p_xyz, fs[0])
                l_coord1, d_coord1 = get_peel_coordinates_of_point_from_face_name(p_xyz, fs[1])
                assert abs(l_coord0 - l_coord1) < 1e-9
                assert abs(d_coord0 - d_coord1) < 1e-9
                l_coord = (l_coord0 + l_coord1)/2
                d_coord = (d_coord0 + d_coord1)/2
            else:
                raise ValueError(f"bad direction {direction}")
        else:
            print(Exception(f"bad number of faces bordering point: {fs}"))

    l_coord = round_off_unwanted_float_precision(l_coord)
    d_coord = round_off_unwanted_float_precision(d_coord)
    return spc, l_coord, d_coord


def round_off_unwanted_float_precision(x):
    # if x has a huge gap in order of magnitude between its second-least-significant digit and its least-significant
    # then the last one is likely float crap
    eps = 1e-9

    if x == 0:
        return x
    if x % eps == 0:
        return x
    
    a = 1
    while x < 1:
        a *= 2
        x *= 2
    
    rem = x % eps
    rem_over_eps = rem / eps
    # print(x, rem_over_eps)
    if rem_over_eps > 1 - 1e-4:
        # we need to ADD a small amount, round x up
        x2 = round(x, 9)
        assert x2 >= x
    elif rem_over_eps < 1e-4:
        x2 = round(x, 9)
        assert x2 <= x
    else:
        # keep the precision
        x2 = x
    return x2/a  # rescale since we multiplied it up


def get_vector_decomposition_coefficients(v, v1, v2):
    # v = a1*v1 + a2*v2, solve for a1 and a2
    x1, y1, z1 = v1
    x2, y2, z2 = v2
    xp, yp, zp = v
    # solve equation for a1 and a2: [[xp] [yp]] = [[x1 x2] [y1 y2]] [[a1] [a2]]
    A = np.array([[x1, x2], [y1, y2]])
    Ainv = np.linalg.inv(A)
    a1, a2 = Ainv @ np.array([xp, yp])
    # verify solution works with z coordinates
    diff = zp - (a1*z1 + a2*z2)
    assert abs(diff) < 1e-9, "z coordinate verification of vector decomposition failed"

    # floats
    if -1e-9 < a1 < 0:
        a1 = 0
    if -1e-9 < a2 < 0:
        a2 = 0

    assert 0 <= a1 <= 1 and 0 <= a2 <= 1, f"vector decomposition should have coefficients between 0 and 1\ngot:\n  {v}\n= {a1} * {v1}\n+ {a2} * {v2}"
    a1 = float(a1)
    a2 = float(a2)
    return a1, a2


def get_peel_coordinates_of_point_from_face_corners(xyz_proj, xyz0, xyz1, xyz2, xyz3):
    # point in question is at xyz_proj AFTER ALREADY BEING PROJECTED onto the plane containing the three vertices of the face it's on
    # face corners are xyz0 at the ancestral starting point, and xyz1/2/3 at the starting points in the 1/2/3 directions from there
    if xyz3 is None:
        d01 = xyz1 - xyz0
        d02 = xyz2 - xyz0
        d0p = xyz_proj - xyz0
        a1, a2 = get_vector_decomposition_coefficients(d0p, d01, d02)

        # 1 direction = L
        # 2 direction = DL
        l_coord = a1 + a2
        d_coord = a2
        assert 0 <= l_coord < 1 and 0 <= d_coord < 1, "left and down coordinates should be between 0 and 1"
    elif xyz1 is None:
        d02 = xyz2 - xyz0
        d03 = xyz3 - xyz0
        d0p = xyz_proj - xyz0
        a2, a3 = get_vector_decomposition_coefficients(d0p, d02, d03)

        # 2 direction = DL
        # 3 direction = D
        l_coord = a2
        d_coord = a2 + a3
        assert 0 <= l_coord < 1 and 0 <= d_coord < 1, "left and down coordinates should be between 0 and 1"
    else:
        raise Exception(f"bad face vertices")
    return l_coord, d_coord


def get_peel_coordinates_of_point_from_face_name(p_xyz, face):
    ax, ay, az, c = fc.get_plane_parameters_of_faces()[face]
    xyz0, xyz1, xyz2, xyz3 = fc.get_face_corner_coordinates_xyz(as_array=True)[face]
    xyz_proj = mcm.project_point_onto_plane(p_xyz, ax, ay, az, c)
    l_coord, d_coord = get_peel_coordinates_of_point_from_face_corners(xyz_proj, xyz0, xyz1, xyz2, xyz3)
    return l_coord, d_coord


def get_vertices_in_common_to_faces(fs):
    s = set(fs[0]) - {"X"}
    for f in fs[1:]:
        s &= set(f)
    return s


def get_point_code_from_peel_coordinates(starting_pc, l_coord, d_coord, max_iterations=15):
    s = starting_pc
    iterations = 0
    while (iterations < max_iterations) and (l_coord > 0 or d_coord > 0):
        # print(f"iteration {iterations}, {l_coord = }, {d_coord = }, {s = }")
        assert 0 <= l_coord < 1 and 0 <= d_coord < 1
        l_coord *= 2
        d_coord *= 2
        l_bit, l_coord = divmod(l_coord, 1)
        d_bit, d_coord = divmod(d_coord, 1)
        assert l_bit in [0,1]
        assert d_bit in [0,1]
        if l_bit == 0 and d_bit == 0:
            c = "0"
        elif l_bit == 0 and d_bit == 1:
            c = "3"
        elif l_bit == 1 and d_bit == 0:
            c = "1"
        elif l_bit == 1 and d_bit == 1:
            c = "2"
        else:
            raise Exception("impossible")
        s += c
        iterations += 1

    # strip trailing zeros which never led to more precision because we hit iteration limit
    while s[-1] == "0":
        s = s[:-1]
    return s


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
    pc = get_point_code_from_peel_coordinates(spc, l, d)
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
    
    spc, l, d = get_peel_coordinates_of_point(p)
    pc = get_point_code_from_peel_coordinates(spc, l, d)
    ll = p.latlondeg(as_array=False)
    fs = fc.get_faces_of_point_by_closest_center(p)
    print(f"latlon {ll}\nis at peel coords L={l}, D={d} from point {spc}\ngot point code {pc} for this point")
    input("check")
    print()
