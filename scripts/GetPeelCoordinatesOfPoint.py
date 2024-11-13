import math
import numpy as np

import icosalattice.IcosahedronMath as icm
import icosalattice.MapCoordinateMath as mcm
import icosalattice.StartingPoints as sp
from icosalattice.UnitSpherePoint import UnitSpherePoint
import icosalattice.Faces as fc


faces = fc.get_face_names()
starting_points, adj = sp.get_starting_points_immutable()
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
    print(p, fs)

    if len(fs) == 1:
        face ,= fs
        ax, ay, az, c = fc.get_plane_parameters_of_faces()[face]
        xyz0, xyz1, xyz2, xyz3 = fc.get_face_corner_coordinates_xyz(as_array=True)[face]
        xyz_proj = mcm.project_point_onto_plane(p_xyz, ax, ay, az, c)
        if xyz3 is None:
            assert fc.get_directionality_of_face(face) == "up"
            d01 = xyz1 - xyz0
            d02 = xyz2 - xyz0
            d0p = xyz_proj - xyz0
            x1, y1, z1 = d01
            x2, y2, z2 = d02
            xp, yp, zp = d0p
            # solve equation for a1 and a2: [[xp] [yp]] = [[x1 x2] [y1 y2]] [[a1] [a2]]
            A = np.array([[x1, x2], [y1, y2]])
            Ainv = np.linalg.inv(A)
            a1, a2 = Ainv @ np.array([xp, yp])
            # verify solution works with z coordinates
            diff = zp - (a1*z1 + a2*z2)
            assert abs(diff) < 1e-9, "z coordinate verification of vector decomposition failed"
            assert 0 <= a1 <= 1 and 0 <= a2 <= 1, "vector decomposition should have coefficients between 0 and 1"
            print(a1, a2)

            # 1 direction = L
            # 2 direction = DL
            l_coord = a1 + a2
            d_coord = a2
            assert 0 <= l_coord <= 1 and 0 <= d_coord <= 1, "left and down coordinates should be between 0 and 1"
        elif xyz1 is None:
            assert fc.get_directionality_of_face(face) == "down"
            d02 = xyz2 - xyz0
            d03 = xyz3 - xyz0
            d0p = xyz_proj - xyz0
            x2, y2, z2 = d02
            x3, y3, z3 = d03
            xp, yp, zp = d0p
            # solve equation for a2 and a3: [[xp] [yp]] = [[x2 x3] [y2 y3]] [[a2] [a3]]
            A = np.array([[x2, x3], [y2, y3]])
            Ainv = np.linalg.inv(A)
            a2, a3 = Ainv @ np.array([xp, yp])
            # verify solution works with z coordinates
            diff = zp - (a2*z2 + a3*z3)
            assert abs(diff) < 1e-9, "z coordinate verification of vector decomposition failed"
            assert 0 <= a2 <= 1 and 0 <= a3 <= 1, "vector decomposition should have coefficients between 0 and 1"
            print(a2, a3)

            # 2 direction = DL
            # 3 direction = D
            l_coord = a2
            d_coord = a2 + a3
            assert 0 <= l_coord <= 1 and 0 <= d_coord <= 1, "left and down coordinates should be between 0 and 1"
        else:
            raise Exception(f"bad face vertices for {face}")
        
        return l_coord, d_coord
    else:
        raise NotImplementedError("multiple faces")

while True:
    # pick a test point
    p = UnitSpherePoint.random()
    print(get_peel_coordinates_of_point(p))