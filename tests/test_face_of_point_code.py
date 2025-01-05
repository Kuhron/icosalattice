import pytest

import math
import numpy as np

import icosalattice.IcosahedronMath as icm
import icosalattice.MapCoordinateMath as mcm
import icosalattice.StartingPoints as sp
import icosalattice.UnitSpherePoint as usp
import icosalattice.Faces as fc
import icosalattice.Edges as ed
from icosalattice.CoordinatesByPlaneGridding import get_xyz_from_point_code_using_corrected_plane_gridding



def test_face_of_point_code():
    faces = fc.get_face_names()
    starting_points, adj = sp.STARTING_POINTS_AND_ADJACENCY
    labels = sp.STARTING_POINT_CODES

    test_cases = {
        "G202202": ["GAEX", "GXEF"],
        "G230031": ["GXEF"],
        "L0222": ["LKJX", "LXJB"],
        "C202113": ["CAKX"],
        "K11011": ["KAIX", "CAKX"],
        "F021": ["FEDX"],
        "F331312": ["FXDB"],
        "I232111": ["IXGH"],
        "H31111": ["HXFB"],
        "C222223": ["CXKL"],
        "E303033": ["EXCD", "FEDX"],
    }
    for pc in sp.STARTING_POINT_CODES:
        test_cases[pc] = sorted([x for x in faces if pc in x])
    for spc, d in sp.STARTING_DIRECTIONAL_DICT.items():
        for x, other_pc in d.items():
            pc = spc + x  # e.g. "C1"
            test_cases[pc] = sorted([x for x in faces if (spc in x and other_pc in x)])

    for pc, faces_expected in test_cases.items():
        faces_expected = sorted(faces_expected)
        fs_from_pc = sorted(fc.get_faces_of_point_code(pc))
        xyz = get_xyz_from_point_code_using_corrected_plane_gridding(pc)
        fs_from_xyz = sorted(fc.get_faces_of_xyz_by_closest_center(xyz))
        assert faces_expected == fs_from_pc == fs_from_xyz, f"expected: {faces_expected} for {pc}\ngot1: {fs_from_pc = }\ngot2: {fs_from_xyz = }"
