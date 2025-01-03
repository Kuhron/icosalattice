import icosalattice.StartingPoints as sp
import icosalattice.MapCoordinateMath as mcm
import icosalattice.PointCodeArithmetic as pca


FACE_NAMES = [
    "CAKX", "CXKL",
    "DCLX", "DXLB",
    "EACX", "EXCD",
    "FEDX", "FXDB",
    "GAEX", "GXEF",
    "HGFX", "HXFB",
    "IAGX", "IXGH",
    "JIHX", "JXHB",
    "KAIX", "KXIJ",
    "LKJX", "LXJB",
]


def get_face_names():
    return FACE_NAMES


def get_face_corner_coordinates_xyz(as_array=False):
    starting_points, adj = sp.STARTING_POINTS_AND_ADJACENCY
    labels = sp.STARTING_POINT_CODES
    label_to_xyz = {label: p.xyz(as_array=as_array) for label, p in zip(labels, starting_points)}
    face_name_to_xyzs = {}
    face_names = get_face_names()
    for face_name in face_names:
        xyzs = [label_to_xyz.get(pc) for pc in face_name]
        face_name_to_xyzs[face_name] = xyzs
    return face_name_to_xyzs


def get_directionality_of_face(face_name):
    # is it an up-pointing triangle or a down-pointing triangle
    if face_name not in FACE_NAMES:
        raise ValueError(f"invalid face name {face_name!r}")
    i = face_name.index("X")
    if i == 1:
        return "down"
    elif i == 3:
        return "up"
    else:
        raise ValueError("bad face name")


def get_plane_parameters_of_faces():
    # ax*x + ay*y + az*z = c
    face_name_to_xyzs = get_face_corner_coordinates_xyz()
    d = {}
    for face_name, xyzs in face_name_to_xyzs.items():
        vertices = [x for x in xyzs if x is not None]
        ax, ay, az, c = mcm.get_plane_containing_three_points_3d(*vertices)
        d[face_name] = (ax, ay, az, c)
    return d


def get_faces_of_point_by_plane_projection(p):
    p_xyz = p.xyz(as_array=True)
    chosen_faces = []
    best_ratio = None
    plane_parameters_by_face = get_plane_parameters_of_faces()
    for face_name in FACE_NAMES:
        ax, ay, az, c = plane_parameters_by_face[face_name]
        ratio = mcm.get_projection_dilation_ratio_of_point_onto_plane(p_xyz, ax, ay, az, c)
        
        if ratio < 0:
            # face is on the other side of the planet, ignore
            continue
        elif ratio > 0:
            if best_ratio is None or abs(ratio - best_ratio) < 1e-9:
                best_ratio = ratio
                chosen_faces.append(face_name)
            elif ratio < best_ratio and abs(ratio-best_ratio) > 1e-9:
                # less dilation of the point's displacement vector means it's closer to this face (TODO not sure I believe this, shouldn't the ratio only be < 1 for one face?)
                best_ratio = ratio
                chosen_faces = [face_name]
            else:
                continue
        else:
            raise ValueError("should not have ratio of zero")

    # p_projected = p_xyz * best_ratio
    # print("by plane projection:", chosen_faces, best_ratio, p_projected)
    return sorted(chosen_faces)


def get_faces_of_point_by_closest_center(p):
    p_xyz = p.xyz(as_array=True)
    return get_faces_of_xyz_by_closest_center(p_xyz)


def get_faces_of_xyz_by_closest_center(xyz):
    face_name_to_xyzs = get_face_corner_coordinates_xyz(as_array=True)
    best_distance = None
    chosen_faces = []
    for face_name, xyzs in face_name_to_xyzs.items():
        xyzs = [x for x in xyzs if x is not None]
        center = sum(xyzs) / 3
        d = mcm.xyz_distance(center, xyz)
        if best_distance is None or (d < best_distance and abs(d-best_distance) > 1e-9):
            # print(f"best distance: {best_distance} -> {d}")
            best_distance = d
            chosen_faces = [face_name]
        elif abs(d - best_distance) < 1e-9:
            chosen_faces.append(face_name)
    # print("by closest center:", chosen_faces, best_distance, p_xyz)
    return sorted(chosen_faces)


def get_vertices_in_common_to_faces(fs):
    s = set(fs[0]) - {"X"}
    for f in fs[1:]:
        s &= set(f)
    return s


def get_faces_in_watershed_of_starting_point(spc):
    fs = [x for x in FACE_NAMES if x.startswith(spc)]
    assert len(fs) == 2, fs
    return fs


def get_faces_of_point_code(pc):
    pc = "".join(x for x in pc if x != "0")  # removing all zeros shouldn't affect the faces it's on
    if len(pc) == 1:
        assert pc in sp.STARTING_POINT_CODES
        return [f for f in FACE_NAMES if pc in f]
    
    spc = pc[0]
    nums = pc[1:]
    spd = sp.STARTING_DIRECTIONAL_DICT
    if len(set(nums)) == 1:
        # only went in one direction, so point is on an edge
        n = nums[0]
        edge = spc + spd[spc][n]
        return [f for f in FACE_NAMES if edge[0] in f and edge[1] in f]
    else:
        # otherwise, the first non-2 direction determines the face (which you can no longer leave once you're on it)
        nums = [x for x in nums if x != "2"]
        x = nums[0]
        if x == "1":
            p0 = spc
            p1 = spd[spc]["1"]
            p2 = spd[spc]["2"]
            p3 = "X"
        elif x == "3":
            p0 = spc
            p1 = "X"
            p2 = spd[spc]["2"]
            p3 = spd[spc]["3"]
        else:
            raise ValueError(f"bad direction: {x}")
        f = p0 + p1 + p2 + p3
        return [f]


def select_point_codes_on_face(pcs, face_name):
    return [pc for pc in pcs if face_name in get_faces_of_point_code(pc)]
