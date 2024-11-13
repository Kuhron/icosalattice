import math
import numpy as np

import icosalattice.IcosahedronMath as icm
import icosalattice.MapCoordinateMath as mcm
import icosalattice.StartingPoints as sp
import icosalattice.UnitSpherePoint as usp
import icosalattice.Faces as fc


faces = fc.get_face_names()
starting_points, adj = sp.get_starting_points_immutable()
labels = sp.STARTING_POINT_CODES
label_to_latlon = {label: p.latlondeg() for label, p in zip(labels, starting_points)}

for face_name, (ax, ay, az, c) in fc.get_plane_parameters_of_faces().items():
    print(f"{face_name}: {ax}*x + {ay}*y + {az}*z = {c}")


# pick a test point
p = usp.UnitSpherePoint.random()
# p = starting_points[2]
f1 = fc.get_faces_of_point_by_plane_projection(p)
f2 = fc.get_faces_of_point_by_closest_center(p)
assert f1 == f2
print(p, f1)

# TODO get coordinates within face/peel
raise NotImplementedError
