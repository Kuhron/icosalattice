import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from icosalattice.CoordinatesOfPointCode import get_xyz_from_point_code, get_latlon_from_point_code, METHOD_NAMES
from icosalattice.GeneratePointCodes import get_all_point_codes_from_ancestor_at_iteration
from icosalattice.Iterations import get_iteration_number_from_point_code


if __name__ == "__main__":
    # pc0 = "C21013"
    # pc0_iterations = get_iteration_number_from_point_code(pc0)
    # pcs = get_all_point_codes_from_ancestor_at_iteration(ancestor_pc=pc0, iterations=pc0_iterations+2)

    pcs = ["C21013", "K22", "I", "D333", "J110102"]

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
