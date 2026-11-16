import random
import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from icosalattice.CoordinatesOfPointCode import get_xyz_from_point_code, get_latlon_from_point_code, METHOD_NAMES
from icosalattice.GeneratePointCodes import get_all_point_codes_from_ancestor_at_iteration
from icosalattice.Iterations import get_iteration_number_from_point_code
from icosalattice.RandomUtil import get_pseudorandom_point_codes
from icosalattice.UnitSpherePoint import UnitSpherePoint
from icosalattice.FindNearestPointCode import get_nearest_point_code_to_latlon


def inspect_output_coordinates():
    # pc0 = "C21013"
    # pc0_iterations = get_iteration_number_from_point_code(pc0)
    # pcs = get_all_point_codes_from_ancestor_at_iteration(ancestor_pc=pc0, iterations=pc0_iterations+2)

    # pcs = ["C21013", "K22", "I", "D333", "J110102"]
    pcs = get_pseudorandom_point_codes(n_pcs=8, seed=1934.85702 + time.time())

    df = pd.DataFrame(index=pcs, columns=[coord+"_"+method for coord in ["x", "y", "z", "lat", "lon"] for method in METHOD_NAMES])
    for pc in pcs:
        for method in METHOD_NAMES:
            x, y, z = get_xyz_from_point_code(pc, method)
            lat, lon = get_latlon_from_point_code(pc, method)
            df.loc[pc, "x_"+method] = x
            df.loc[pc, "y_"+method] = y
            df.loc[pc, "z_"+method] = z
            df.loc[pc, "lat_"+method] = lat
            df.loc[pc, "lon_"+method] = lon
    print(df.T)


def compare_speeds():
    # speed comparison
    print("pc -> latlon")
    pcs = get_pseudorandom_point_codes(n_pcs=5000, seed=4444.77073918 + time.time())
    for method in METHOD_NAMES:
        print(f"running method {method}")
        t0 = time.time()
        for pc in pcs:
            latlon = get_latlon_from_point_code(pc, method=method)
        dt = time.time() - t0
        print(f"method {method} took {dt:.6f} seconds")
    print()

    print("latlon -> pc")
    latlons = []
    for i in range(1000):
        usp = UnitSpherePoint.get_random_unit_sphere_point()
        latlon = usp.latlondeg()
        latlons.append(latlon)

    for method in METHOD_NAMES:
        print(f"running method {method}")
        t0 = time.time()
        for latlon in latlons:
            nearest_candidate_pc, distance_normalized, distance_in_units = get_nearest_point_code_to_latlon(latlon, maximum_distance=1, planet_radius=6371, pc_to_xyz_method=method)
        dt = time.time() - t0
        print(f"method {method} took {dt:.6f} seconds")



if __name__ == "__main__":
    inspect_output_coordinates()
    print()
    compare_speeds()

    # ebs1 is closest to what I'd consider the "ground truth", since it's the direct implementation of how this whole coordinate system was conceived (iterative bisecting of edges starting from icosahedron)
    # cpg1 is quite similar to ebs1 and much faster to compute (at least for pc -> coords, since they're all slow with the coords -> pc search)
    # upg1 is the fastest to compute but only like a factor of 2 faster than cpg1, and it has significantly more distortion so it's not worth the computational gain
    # rta1 is just doing its own weird thing, don't use it

    # so cpg1 is the clear choice
    # - should be the default method until a better one is found (if I can ever figure out the Jacobian stuff to get the ideal adjustment function)
