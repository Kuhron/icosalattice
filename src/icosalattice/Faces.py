import icosalattice.StartingPoints as sp
import icosalattice.MapCoordinateMath as mcm
import warnings


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
    face_name_to_xyzs = get_face_corner_coordinates_xyz(as_array=True)
    best_distance = None
    chosen_faces = []
    for face_name, xyzs in face_name_to_xyzs.items():
        xyzs = [x for x in xyzs if x is not None]
        center = sum(xyzs) / 3
        d = mcm.xyz_distance(center, p_xyz)
        if best_distance is None or (d < best_distance and abs(d-best_distance) > 1e-9):
            # print(f"best distance: {best_distance} -> {d}")
            best_distance = d
            chosen_faces = [face_name]
        elif abs(d - best_distance) < 1e-9:
            chosen_faces.append(face_name)
    # print("by closest center:", chosen_faces, best_distance, p_xyz)
    return sorted(chosen_faces)
