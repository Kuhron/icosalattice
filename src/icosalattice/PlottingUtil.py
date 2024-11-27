# for more basic plotting functions
# more specific things about certain kinds of plot (e.g. adjacency directions, point codes on half-peels, etc.) should be in specialized modules


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CloughTocher2DInterpolator


def plot_interpolated_data(xys, zs, xs_of_grid, ys_of_grid, xlim=None, ylim=None, title=None, show=True):
    # interpolate from the data (zs) located at the points (xys) onto the grid given by x_range, y_range, n_xs, n_ys
    X, Y = np.meshgrid(xs_of_grid, ys_of_grid)
    interp = CloughTocher2DInterpolator(xys, zs)
    Z = interp(X, Y)
    plt.pcolormesh(X, Y, Z, shading='auto')
    plt.colorbar()
    plt.axis("equal")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.title(title)
    if show:
        plt.show()
