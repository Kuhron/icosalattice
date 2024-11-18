import pytest

import numpy as np

import icosalattice.IcosahedronMath as icm
import icosalattice.MapCoordinateMath as mcm
import icosalattice.StartingPoints as sp
import icosalattice.Faces as fc


def test_face_areas():
    starting_points, adj = sp.STARTING_POINTS_AND_ADJACENCY
    labels = sp.STARTING_POINT_CODES
    
    face_name_to_xyzs = fc.get_face_corner_coordinates_xyz()

    areas = set()
    for face_name, xyzs in face_name_to_xyzs.items():
        vertices = [x for x in xyzs if x is not None]
        area = mcm.area_of_triangle_from_vertices_3d(*vertices)
        # print(face_name, area)
        areas.add(area)
    assert max(areas) - min(areas) < 1e-9, "unequal face areas"
    area = np.mean(list(areas))

    # have another test for edge lengths being correct, so just get one of them here
    edge_length = mcm.xyz_distance(starting_points[labels.index("A")].xyz(), starting_points[labels.index("C")].xyz())
    expected_surface_area = 5 * (3**0.5) * edge_length**2  # https://en.wikipedia.org/wiki/Regular_icosahedron#Properties
    got_surface_area = 20 * area
    assert abs(expected_surface_area - got_surface_area) < 1e-9, "bad surface area"
