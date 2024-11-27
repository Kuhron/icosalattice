import numpy as np
import matplotlib.pyplot as plt

import icosalattice.IcosahedronMath as icm
from icosalattice.PlottingUtil import plot_interpolated_data
from icosalattice.CoordinatesOfPointCode import get_latlon_from_point_code



def plot_variable_interpolated_from_dict(pc_to_val, dots_per_degree, title=None, show=True):
    pcs, vals = zip(*pc_to_val.items())
    latlons = [get_latlon_from_point_code(pc) for pc in pcs]
    xys = convert_latlons_to_xys_for_equirectangular_projection(latlons)
    xs_of_grid = np.linspace(-180, 180, 360*dots_per_degree)
    ys_of_grid = np.linspace(-90, 90, 180*dots_per_degree)
    xlim = (-180, 180)
    ylim = (-90, 90)
    plot_interpolated_data(xys, vals, xs_of_grid, ys_of_grid, xlim=xlim, ylim=ylim, title=title, show=show)


def convert_latlons_to_xys_for_equirectangular_projection(latlons):
    lats, lons = zip(*latlons)
    ys = lats
    xs = lons
    xys = list(zip(xs, ys))
    return xys
