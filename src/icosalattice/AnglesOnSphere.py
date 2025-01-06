# functions relating to measuring angles on the sphere surface

import math
import numpy as np

from icosalattice.DistancesOnSphere import distance_great_circle


# TODO will probably be good to implement generic methods for solving spherical triangles
# TODO implement measurement of area occupied by a set of points that approximate a 2D region, e.g. a country; use sum of a bunch of triangles for this?

# https://en.wikipedia.org/wiki/Spherical_trigonometry


def measure_angles_on_sphere(xyz_A, xyz_B, xyz_C):
    # angle ABC with vertex at B; what angle do the paths AB and BC make when they meet at B?
    # assumes radius = 1
    a = distance_great_circle(xyz_B, xyz_C)
    b = distance_great_circle(xyz_A, xyz_C)
    c = distance_great_circle(xyz_A, xyz_B)
    cos_a = np.cos(a)
    cos_b = np.cos(b)
    cos_c = np.cos(c)
    sin_a = np.sin(a)
    sin_b = np.sin(b)
    sin_c = np.sin(c)
    cos_A = (cos_a - cos_b*cos_c)/(sin_b*sin_c)
    cos_B = (cos_b - cos_c*cos_a)/(sin_c*sin_a)
    cos_C = (cos_c - cos_a*cos_b)/(sin_a*sin_b)
    return [np.arccos(cos_A), np.arccos(cos_B), np.arccos(cos_C)]


def get_area_of_triangle_on_sphere(xyz1, xyz2, xyz3):
    # Girard's theorem
    # assumes radius = 1
    angles = measure_angles_on_sphere(xyz1, xyz2, xyz3)
    return sum(angles) - np.pi

