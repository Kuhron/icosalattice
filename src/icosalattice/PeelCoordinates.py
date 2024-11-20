from icosalattice.UnitSpherePoint import UnitSpherePoint
import icosalattice.MathUtil as mu
import icosalattice.StartingPoints as sp
import icosalattice.Faces as fc
import icosalattice.Edges as ed
import icosalattice.MapCoordinateMath as mcm
import numpy as np



def get_peel_coordinates_of_point(p: UnitSpherePoint):
    p_xyz = p.xyz(as_array=True)
    return get_peel_coordinates_from_xyz(p_xyz)


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
    spc = pc[0]
    l_coord = 0
    d_coord = 0
    c_to_bits = {
        "0": (0, 0), "1": (1, 0),
        "3": (0, 1), "2": (1, 1),
    }
    for i, c in enumerate(pc[1:]):
        denom = 2**(i+1)
        l_bit, d_bit = c_to_bits[c]
        l_coord += l_bit/denom
        d_coord += d_bit/denom
    return spc, l_coord, d_coord


def get_xyz_from_peel_coordinates(starting_pc, l_coord, d_coord):
    if starting_pc in ["A", "B"]:
        assert l_coord == 0 and d_coord == 0, f"cannot have steps from a pole, but got spc={starting_pc}, l={l_coord}, d={d_coord}"

    assert 0 <= l_coord < 1 and 0 <= d_coord < 1, f"must have peel coords in interval [0, 1), but got l={l_coord}, d={d_coord}"

    if l_coord == 0 and d_coord == 0:
        return sp.STARTING_POINTS[sp.STARTING_POINT_CODES.index(starting_pc)].xyz()
    
    fs = [x for x in fc.FACE_NAMES if x.startswith(starting_pc)]
    assert len(fs) == 2, fs
    face_name_to_xyzs = fc.get_face_corner_coordinates_xyz(as_array=True)

    if l_coord > d_coord:
        # on upward-pointing face of this half-peel
        f ,= [x for x in fs if fc.get_directionality_of_face(x) == "up"]
        vertices_xyz = face_name_to_xyzs[f]
        return get_xyz_from_peel_coordinates_on_upward_face(vertices_xyz, l_coord, d_coord)
    elif d_coord > l_coord:
        # on downward-pointing face of this half-peel
        f ,= [x for x in fs if fc.get_directionality_of_face(x) == "down"]
        vertices_xyz = face_name_to_xyzs[f]
        return get_xyz_from_peel_coordinates_on_downward_face(vertices_xyz, l_coord, d_coord)
    else:
        # test both faces
        vertices_xyz1 = face_name_to_xyzs[fs[0]]
        vertices_xyz2 = face_name_to_xyzs[fs[1]]
        xyz1 = get_xyz_from_peel_coordinates_on_upward_face(vertices_xyz1, l_coord, d_coord)
        xyz2 = get_xyz_from_peel_coordinates_on_downward_face(vertices_xyz2, l_coord, d_coord)
        assert (abs(xyz1 - xyz2) < 1e-9).all(), f"{xyz1} != {xyz2}"
        return (xyz1 + xyz2)/2


def get_xyz_from_peel_coordinates_on_upward_face(vertices_xyz, l_coord, d_coord):
    assert l_coord >= d_coord, f"upward face should have l >= d, but got l = {l_coord}, d = {d_coord}"
    xyz0, xyz1, xyz2, xyz3 = vertices_xyz
    assert xyz3 is None, f"vertices do not respect upward-pointing face polarity: {vertices_xyz}"

    l_vector = xyz1 - xyz0
    dl_vector = xyz2 - xyz0

    dl_coefficient = d_coord
    l_coefficient = l_coord - d_coord  # whatever L component is leftover after removing common movement to both D and L
    assert 0 <= dl_coefficient < 1
    assert 0 <= l_coefficient < 1

    dxyz_on_plane = l_coefficient * l_vector + dl_coefficient * dl_vector
    xyz_on_plane = xyz0 + dxyz_on_plane
    xyz_on_sphere = xyz_on_plane / np.linalg.norm(xyz_on_plane)
    return xyz_on_sphere


def get_xyz_from_peel_coordinates_on_downward_face(vertices_xyz, l_coord, d_coord):
    assert d_coord >= l_coord, f"downward face should have d >= l, but got l = {l_coord}, d = {d_coord}"
    xyz0, xyz1, xyz2, xyz3 = vertices_xyz
    assert xyz1 is None, f"vertices do not respect downward-pointing face polarity: {vertices_xyz}"

    dl_vector = xyz2 - xyz0
    d_vector = xyz3 - xyz0

    dl_coefficient = l_coord
    d_coefficient = d_coord - l_coord  # whatever D component is leftover after removing common movement to both D and L
    assert 0 <= dl_coefficient < 1
    assert 0 <= d_coefficient < 1

    dxyz_on_plane = dl_coefficient * dl_vector + d_coefficient * d_vector
    xyz_on_plane = xyz0 + dxyz_on_plane
    xyz_on_sphere = xyz_on_plane / np.linalg.norm(xyz_on_plane)
    return xyz_on_sphere


def get_peel_coordinates_from_xyz(xyz):
    # peel coordinates are within the shape formed by combining this face and its up-facing or down/facing counterpart
    # and shearing those such that it's a square with the 0 point (C,D,E,F,G,H,I,J,K,L) on the upper right
    # the 1 direction going left, the 2 direction going down-left, and the 3 direction going down
    # so each peel has a north and a south square, each the parental watershed of one of the 10 mid-latitude starting points
    # within this square, measure how far left we go and how far down we go (each in interval from 0 to 1)
    # from this we should be able to get the point code directly, I hope?
    fs = fc.get_faces_of_xyz_by_closest_center(xyz)

    if len(fs) == 1:
        face ,= fs
        spc = face[0]
        l_coord, d_coord = get_peel_coordinates_of_point_from_face_name(xyz, face)
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
                s = fc.get_vertices_in_common_to_faces(fs)
                assert len(s) == 1, s
                spc ,= s
            l_coord = 0
            d_coord = 0
        elif len(fs) == 2:
            s = fc.get_vertices_in_common_to_faces(fs)
            spc, direction = ed.get_ancestor_starting_point_and_direction_of_edge(s)
            # print(f"{spc = }, {direction = }")

            if direction == "1":
                # get peel coordinates on up-pointing face radiating from ancestor point
                face ,= [x for x in fs if x[0] == spc]
                l_coord, d_coord = get_peel_coordinates_of_point_from_face_name(xyz, face)
            elif direction == "3":
                # edge point in the 2 direction, so it borders both the up-pointing and down-pointing faces
                # verify that getting coords from both faces give same result
                face ,= [x for x in fs if x[0] == spc]
                l_coord, d_coord = get_peel_coordinates_of_point_from_face_name(xyz, face)
            elif direction == "2":
                # get peel coordinates on down-pointing face radiating from ancestor point
                l_coord0, d_coord0 = get_peel_coordinates_of_point_from_face_name(xyz, fs[0])
                l_coord1, d_coord1 = get_peel_coordinates_of_point_from_face_name(xyz, fs[1])
                assert abs(l_coord0 - l_coord1) < 1e-9
                assert abs(d_coord0 - d_coord1) < 1e-9
                l_coord = (l_coord0 + l_coord1)/2
                d_coord = (d_coord0 + d_coord1)/2
            else:
                raise ValueError(f"bad direction {direction}")
        else:
            raise Exception(f"bad number of faces bordering point: {fs}")

    l_coord = mu.round_off_unwanted_float_precision(l_coord)
    d_coord = mu.round_off_unwanted_float_precision(d_coord)
    return spc, l_coord, d_coord


def get_xyz_from_point_code_using_peel_coordinates(pc):
    spc, l, d = get_peel_coordinates_from_point_code(pc)
    xyz = get_xyz_from_peel_coordinates(spc, l, d)
    return xyz


def get_point_code_from_xyz_using_peel_coordinates(xyz):
    spc, l, d = get_peel_coordinates_from_xyz(xyz)
    pc = get_point_code_from_peel_coordinates(spc, l, d)
    return pc



def get_peel_coordinates_of_point_from_face_corners(xyz_proj, xyz0, xyz1, xyz2, xyz3):
    # point in question is at xyz_proj AFTER ALREADY BEING PROJECTED onto the plane containing the three vertices of the face it's on
    # face corners are xyz0 at the ancestral starting point, and xyz1/2/3 at the starting points in the 1/2/3 directions from there
    if xyz3 is None:
        d01 = xyz1 - xyz0
        d02 = xyz2 - xyz0
        d0p = xyz_proj - xyz0
        a1, a2 = mu.get_vector_decomposition_coefficients(d0p, d01, d02)

        # 1 direction = L
        # 2 direction = DL
        l_coord = a1 + a2
        d_coord = a2
        assert 0 <= l_coord < 1 and 0 <= d_coord < 1, "left and down coordinates should be between 0 and 1"
    elif xyz1 is None:
        d02 = xyz2 - xyz0
        d03 = xyz3 - xyz0
        d0p = xyz_proj - xyz0
        a2, a3 = mu.get_vector_decomposition_coefficients(d0p, d02, d03)

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

