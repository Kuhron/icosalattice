# peel coordinates are within the shape formed by combining this face and its up-facing or down/facing counterpart
# and shearing those such that it's a square with the 0 point (C,D,E,F,G,H,I,J,K,L) on the upper right
# the 1 direction going left, the 2 direction going down-left, and the 3 direction going down
# so each peel has a north and a south square, each the parental watershed of one of the 10 mid-latitude starting points
# within this square, measure how far left we go and how far down we go (each in interval from 0 to 1)
# from this we can get the point code directly

# "raw" peel coordinates means those encoded directly in the point code (which IS reflective of their position on the sphere surface, but not face plane)
# - e.g. C1 is treated as halfway from C to A on the face plane, C11 is treated as halfway from C1 to A on the face plane, etc.
# - as if there is no distortion when projecting between the face plane and the sphere surface
# "adjusted" peel coordinates means those encoding where the point from the sphere surface actually lands on the face plane when projected
# - adjusting moves raw peel coordinates closer to the center, with a few exceptions such as the corners and center of the face, which don't move
# - because points nearer the center of the face plane get "bowed out" / spread out more when projected outward to the sphere surface


import numpy as np

from icosalattice.UnitSpherePoint import UnitSpherePoint
import icosalattice.MathUtil as mu
import icosalattice.StartingPoints as sp
import icosalattice.Faces as fc
import icosalattice.Edges as ed
import icosalattice.MapCoordinateMath as mcm
import icosalattice.CoordinatesByAncestry as anc
import icosalattice.TriangularPeelCoordinates as tri
import icosalattice.PointCodeArithmetic as pca


# potential optimizations, if needed:
# - for the vector composition stuff for getting peel coordinates from xyz and vice versa
# - - could be good to keep all the relevant displacement vectors as constants in one of the modules (maybe Edges)
# - - and ideally have their values analytically determined and specified like how MID_LAT_DEG is known from trig
# - - then can just look them up without recalculating them


def get_adjusted_peel_coordinates_from_raw_peel_coordinates(sld):
    spc, l_raw, d_raw = sld
    l_adjusted, d_adjusted = tri.adjust_ld_using_lp_transformation_in_triangle_coordinates(l_raw, d_raw)
    return spc, l_adjusted, d_adjusted


def get_raw_peel_coordinates_from_adjusted_peel_coordinates(sld):
    spc, l_adjusted, d_adjusted = sld
    l_raw, d_raw = tri.deadjust_ld_using_lp_transformation_in_triangle_coordinates(l_adjusted, d_adjusted)
    return spc, l_raw, d_raw


def get_adjusted_peel_coordinates_of_point(p: UnitSpherePoint):
    p_xyz = p.xyz(as_array=True)
    return get_adjusted_peel_coordinates_from_xyz(p_xyz)


def get_point_code_from_raw_peel_coordinates(sld, max_iterations=15, allow_clipping=False):
    # take the point code at face value as encoding peel coordinates,
    # ignoring any questions about distortion or where exactly the l and d coordinates are located
    starting_pc, l_coord, d_coord = sld
    assert starting_pc in sp.STARTING_POINT_CODES, f"{s} is not a valid starting point code"

    # get bits for l and d first
    iterations = 0
    l_bits = []
    d_bits = []
    while l_coord > 0 or d_coord > 0:
        # print(f"iteration {iterations}, {l_coord = }, {d_coord = }, {s = }")
        if iterations >= max_iterations:
            if allow_clipping:
                print(f"max iterations reached for getting point code from peel coordinates; clipping")
                break
            else:
                raise RuntimeError(f"failed to finish getting point code from peel coordinates without clipping\n{l_coord = }, {d_coord = }, current pc = {s}")
        
        assert 0 <= l_coord < 1 and 0 <= d_coord < 1
        l_coord *= 2
        d_coord *= 2

        l_bit, l_coord = divmod(l_coord, 1)
        d_bit, d_coord = divmod(d_coord, 1)
        assert l_bit in [0,1]
        assert d_bit in [0,1]
        l_bits.append(l_bit)
        d_bits.append(d_bit)

        iterations += 1

    assert len(l_bits) <= max_iterations
    assert len(d_bits) <= max_iterations

    # round up based on the remaining amount
    l_bits = round_bit_array(l_bits, rem=l_coord)
    d_bits = round_bit_array(d_bits, rem=d_coord)

    s = get_point_code_from_bit_arrays(starting_pc, l_bits, d_bits, strip_trailing_zeros=True)
    assert len(s) <= max_iterations + 1, f"{max_iterations = }, but got {s = }"
    return s


def get_point_code_from_bit_arrays(spc, l_bits, d_bits, strip_trailing_zeros=True):
    s = spc
    for l_bit, d_bit in zip(l_bits, d_bits):
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
    if strip_trailing_zeros:
        s = pca.strip_trailing_zeros(s)
    return s


def round_bit_array(bits, rem):
    # rem is what remains that we would have used for the next order of magnitude down
    assert 0 <= rem < 1, f"{rem = }, should have 0 <= rem < 1"
    assert all(x in [0, 1] for x in bits), f"all bits should be 0 or 1, got these values: {sorted(set(bits))}"
    if rem < 0.5:
        return bits
    else:
        new_bits = [x for x in bits]
        i = len(bits) - 1
        while bits[i] == 1:
            new_bits[i] = 0
            i -= 1
            if i < 0:
                raise RuntimeError(f"overflow when rounding up bit array: {bits}")
        new_bits[i] = 1
        return new_bits


def get_raw_peel_coordinates_from_point_code(pc):
    # take the point code at face value as encoding peel coordinates,
    # ignoring any questions about distortion or where exactly the l and d coordinates are located
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


def get_xyz_from_adjusted_peel_coordinates(sld, as_array=True):
    # take l_coord and d_coord at face value as linear proportions of the way through the half-peel
    # project the resulting point from the face plane to the sphere surface
    spc, l, d = sld
    if spc in ["A", "B"]:
        assert l == 0 and d == 0, f"cannot have steps from a pole, but got spc={spc}, l={l}, d={d}"

    assert 0 <= l < 1 and 0 <= d < 1, f"must have peel coords in interval [0, 1), but got l={l}, d={d}"

    if l == 0 and d == 0:
        return sp.STARTING_POINTS[sp.STARTING_POINT_CODES.index(spc)].xyz(as_array=as_array)
    
    fs = fc.get_faces_in_watershed_of_starting_point(spc)
    face_name_to_xyzs = fc.get_face_corner_coordinates_xyz(as_array=True)

    if l > d:
        # on upward-pointing face of this half-peel
        f ,= [x for x in fs if fc.get_directionality_of_face(x) == "up"]
        vertices_xyz = face_name_to_xyzs[f]
        return get_xyz_from_adjusted_peel_coordinates_on_upward_face(vertices_xyz, l, d, as_array=as_array)
    elif d > l:
        # on downward-pointing face of this half-peel
        f ,= [x for x in fs if fc.get_directionality_of_face(x) == "down"]
        vertices_xyz = face_name_to_xyzs[f]
        return get_xyz_from_adjusted_peel_coordinates_on_downward_face(vertices_xyz, l, d, as_array=as_array)
    else:
        # test both faces
        vertices_xyz1 = face_name_to_xyzs[fs[0]]
        vertices_xyz2 = face_name_to_xyzs[fs[1]]
        xyz1 = get_xyz_from_adjusted_peel_coordinates_on_upward_face(vertices_xyz1, l, d, as_array=True)
        xyz2 = get_xyz_from_adjusted_peel_coordinates_on_downward_face(vertices_xyz2, l, d, as_array=True)
        assert (abs(xyz1 - xyz2) < 1e-9).all(), f"{xyz1} != {xyz2}"
        res = (xyz1 + xyz2)/2
        if not as_array:
            x, y, z = res
            return (float(x), float(y), float(z))
        return res


def get_xyz_from_adjusted_peel_coordinates_on_upward_face(vertices_xyz, l, d, as_array=True):
    assert l >= d, f"upward face should have l >= d, but got l = {l}, d = {d}"
    xyz0, xyz1, xyz2, xyz3 = vertices_xyz
    assert xyz3 is None, f"vertices do not respect upward-pointing face polarity: {vertices_xyz}"

    l_vector = xyz1 - xyz0
    dl_vector = xyz2 - xyz0

    dl_coefficient = d
    l_coefficient = l - d  # whatever L component is leftover after removing common movement to both D and L
    assert 0 <= dl_coefficient < 1
    assert 0 <= l_coefficient < 1

    dxyz_on_plane = l_coefficient * l_vector + dl_coefficient * dl_vector
    xyz_on_plane = xyz0 + dxyz_on_plane
    xyz_on_sphere = xyz_on_plane / np.linalg.norm(xyz_on_plane)

    if as_array:
        return xyz_on_sphere
    else:
        x,y,z = xyz_on_sphere
        return (float(x), float(y), float(z))


def get_xyz_from_adjusted_peel_coordinates_on_downward_face(vertices_xyz, l, d, as_array=True):
    assert d >= l, f"downward face should have d >= l, but got l = {l}, d = {d}"
    xyz0, xyz1, xyz2, xyz3 = vertices_xyz
    assert xyz1 is None, f"vertices do not respect downward-pointing face polarity: {vertices_xyz}"

    dl_vector = xyz2 - xyz0
    d_vector = xyz3 - xyz0

    dl_coefficient = l
    d_coefficient = d - l  # whatever D component is leftover after removing common movement to both D and L
    assert 0 <= dl_coefficient < 1
    assert 0 <= d_coefficient < 1

    dxyz_on_plane = dl_coefficient * dl_vector + d_coefficient * d_vector
    xyz_on_plane = xyz0 + dxyz_on_plane
    xyz_on_sphere = xyz_on_plane / np.linalg.norm(xyz_on_plane)
    
    if as_array:
        return xyz_on_sphere
    else:
        x,y,z = xyz_on_sphere
        return (float(x), float(y), float(z))


def get_adjusted_peel_coordinates_from_xyz(xyz):
    fs = fc.get_faces_of_xyz_by_closest_center(xyz)

    if len(fs) == 1:
        face ,= fs
        spc = face[0]
        l_coord, d_coord = get_adjusted_ld_of_point_on_face(xyz, face)
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
                l_coord, d_coord = get_adjusted_ld_of_point_on_face(xyz, face)
            elif direction == "3":
                # get peel coordinates on down-pointing face radiating from ancestor point
                face ,= [x for x in fs if x[0] == spc]
                l_coord, d_coord = get_adjusted_ld_of_point_on_face(xyz, face)
            elif direction == "2":
                # edge point in the 2 direction, so it borders both the up-pointing and down-pointing faces
                # verify that getting coords from both faces give same result
                l_coord0, d_coord0 = get_adjusted_ld_of_point_on_face(xyz, fs[0])
                l_coord1, d_coord1 = get_adjusted_ld_of_point_on_face(xyz, fs[1])
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


def get_xyz_from_point_code_using_peel_coordinates(pc, as_array=True):
    sld_raw = get_raw_peel_coordinates_from_point_code(pc)
    sld_adjusted = get_adjusted_peel_coordinates_from_raw_peel_coordinates(sld_raw)
    xyz = get_xyz_from_adjusted_peel_coordinates(sld_adjusted, as_array=as_array)
    return xyz


def get_point_code_from_xyz_using_peel_coordinates(xyz):
    sld_adjusted = get_adjusted_peel_coordinates_from_xyz(xyz)
    sld_raw = get_raw_peel_coordinates_from_adjusted_peel_coordinates(sld_adjusted)
    pc = get_point_code_from_raw_peel_coordinates(sld_raw)
    return pc


def get_latlon_from_point_code_using_peel_coordinates(pc, deg=True, as_array=True):
    xyz = get_xyz_from_point_code_using_peel_coordinates(pc)
    return mcm.unit_vector_cartesian_to_latlon(*xyz, deg=deg, as_array=as_array)


def get_point_code_from_latlon_using_peel_coordinates(latlon, deg=True, as_array=True):
    xyz = mcm.unit_vector_latlon_to_cartesian(*latlon, deg=deg, as_array=as_array)
    return get_point_code_from_xyz_using_peel_coordinates(xyz)


def get_adjusted_ld_of_point_from_face_corners(xyz_proj, xyz0, xyz1, xyz2, xyz3, allow_one=False):
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
        validate_l_and_d_coordinates(l_coord, d_coord, allow_one=allow_one)
    elif xyz1 is None:
        d02 = xyz2 - xyz0
        d03 = xyz3 - xyz0
        d0p = xyz_proj - xyz0
        a2, a3 = mu.get_vector_decomposition_coefficients(d0p, d02, d03)

        # 2 direction = DL
        # 3 direction = D
        l_coord = a2
        d_coord = a2 + a3
        validate_l_and_d_coordinates(l_coord, d_coord, allow_one=allow_one)
    else:
        raise Exception(f"bad face vertices")
    return l_coord, d_coord


def validate_l_and_d_coordinates(l_coord, d_coord, allow_one=False):
    l_geq_0 = abs(l_coord) < 1e-9 or l_coord > 0
    d_geq_0 = abs(d_coord) < 1e-9 or d_coord > 0
    l_lt_1 = l_coord < 1
    d_lt_1 = d_coord < 1
    l_eq_1 = abs(1 - l_coord) < 1e-9
    d_eq_1 = abs(1 - d_coord) < 1e-9
    if allow_one:
        l_satisfied = l_geq_0 and (l_lt_1 or l_eq_1)
        d_satisfied = d_geq_0 and (d_lt_1 or d_eq_1)
        assert l_satisfied and d_satisfied, f"left and down coordinates should be between 0 and 1 (here, allowing equal to one)\nbut got (l, d) = {(l_coord, d_coord)}"
    else:
        l_satisfied = l_geq_0 and l_lt_1
        d_satisfied = d_geq_0 and d_lt_1
        assert l_satisfied and d_satisfied, f"left and down coordinates should be between 0 and 1 (here, excluding equal to one)\nbut got (l, d) = {(l_coord, d_coord)}"


def get_adjusted_ld_of_point_on_face(p_xyz, face_name, allow_one=False):
    ax, ay, az, c = fc.get_plane_parameters_of_faces()[face_name]
    xyz0, xyz1, xyz2, xyz3 = fc.get_face_corner_coordinates_xyz(as_array=True)[face_name]
    xyz_proj = mcm.project_point_onto_plane(p_xyz, ax, ay, az, c)
    l, d = get_adjusted_ld_of_point_from_face_corners(xyz_proj, xyz0, xyz1, xyz2, xyz3, allow_one=allow_one)
    return l, d


def get_adjusted_peel_coordinates_of_point_codes_on_face(pcs, face_name, func_pc_to_xyz):
    ls = []
    ds = []
    for pc in pcs:
        xyz = func_pc_to_xyz(pc, as_array=True)
        try:
            l, d = get_adjusted_ld_of_point_on_face(xyz, face_name, allow_one=True)
        except mu.InvalidVectorDecompositionException:
            l, d = None, None
        ls.append(l)
        ds.append(d)
    
    return ls, ds
