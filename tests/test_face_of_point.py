import pytest

import math
import numpy as np

import icosalattice.IcosahedronMath as icm
import icosalattice.MapCoordinateMath as mcm
import icosalattice.StartingPoints as sp
import icosalattice.UnitSpherePoint as usp
import icosalattice.Faces as fc
import icosalattice.Edges as ed



def test_face_of_point():
    faces = fc.get_face_names()
    starting_points, adj = sp.STARTING_POINTS_AND_ADJACENCY
    labels = sp.STARTING_POINT_CODES
    # label_to_latlon = {label: p.latlondeg() for label, p in zip(labels, starting_points)}

    test_cases = {
        usp.UnitSpherePoint.from_latlondeg(1, 104): ["GXEF"],
        usp.UnitSpherePoint.from_latlondeg(-12, -77): ["LKJX"],
        usp.UnitSpherePoint.from_latlondeg(64, -22): ["CAKX"],
        usp.UnitSpherePoint.from_latlondeg(47, -122): ["KAIX"],
        usp.UnitSpherePoint.from_latlondeg(-19, 48): ["FEDX"],
        usp.UnitSpherePoint.from_latlondeg(-49, 70): ["FXDB"],
        usp.UnitSpherePoint.from_latlondeg(21, -158): ["IXGH"],
        usp.UnitSpherePoint.from_latlondeg(-43, 147): ["HXFB"],
        usp.UnitSpherePoint.from_latlondeg(15, -17): ["CXKL"],
    }
    for i, pc in enumerate(sp.STARTING_POINT_CODES):
        test_cases[starting_points[i]] = sorted([x for x in faces if pc in x])
    edge_midpoints = ed.get_edge_midpoints()
    for label, p in edge_midpoints.items():
        test_cases[p] = sorted([x for x in faces if (label[0] in x and label[1] in x)])

    for p, faces_expected in test_cases.items():
        f1 = fc.get_faces_of_point_by_plane_projection(p)
        f2 = fc.get_faces_of_point_by_closest_center(p)
        assert faces_expected == f1 == f2, f"expected: {faces_expected} for {p}\ngot1: {f1}\ngot2: {f2}"
