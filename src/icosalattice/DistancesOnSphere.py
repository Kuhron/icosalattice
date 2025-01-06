# functions relating to measuring distance on sphere surface

import math
import numpy as np



def convert_distance_3d_to_great_circle(d0, r=1):
    theta = 2 * np.arcsin(d0 / (2*r))
    d_gc = r * theta
    # assert (0 <= d_gc).all(), f"bad great circle distance {d_gc} from d0={d0}, r={r}"
    # assert (d_gc <= np.pi * r).all(), f"bad great circle distance {d_gc} from d0={d0}, r={r}"
    # assert ((d_gc > d0) | (abs(d_gc - d0) < 1e-9)).all(), f"shortest distance should be a straight line, but got great-circle {d_gc} from Euclidean {d0}"
    # print(f"d0 = {d0}, r = {r} -> great circle distance {d_gc}")
    return d_gc
    # return (np.vectorize(lambda d: UnitSpherePoint.convert_distance_3d_to_great_circle_single_value(d, radius=radius)))(d0)


def convert_distance_great_circle_to_3d(d_gc, r=1):
    theta = d_gc / r
    d0 = 2 * r * np.sin(theta)
    assert 0 <= d0 <= 2*r, f"bad 3d distance {d0} from d_gc={d_gc}, r={r}"
    assert d0 <= d_gc, "shortest distance should be a straight line"
    return d0


def distance_3d(xyz1, xyz2):
    return np.linalg.norm(xyz1 - xyz2)


def distance_great_circle(xyz1, xyz2, r=1):
    d0 = distance_3d(xyz1, xyz2)
    return convert_distance_3d_to_great_circle(d0, r)
