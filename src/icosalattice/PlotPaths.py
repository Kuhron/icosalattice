# looking at paths of successive points to see how the direction changes, how far apart each pair of points is, etc.

import numpy as np
import matplotlib.pyplot as plt


def plot_distances_and_angles_2d(distances, angles):
    plt.subplot(2, 1, 1)
    plt.title("in 2D")
    plt.plot(distances)
    plt.ylabel("distance")
    plt.subplot(2, 1, 2)
    plt.plot(angles * 180/np.pi)
    plt.ylabel("deg")
    plt.show()


def plot_distances_and_angles_3d(distances, angles_xy, angles_xz, angles_yz):
    plt.subplot(2, 1, 1)
    plt.title("in 3D")
    plt.plot(distances)
    plt.ylabel("distance")
    plt.subplot(2, 1, 2)
    plt.plot(angles_xy * 180/np.pi, label="xy")
    plt.plot(angles_xz * 180/np.pi, label="xz")
    plt.plot(angles_yz * 180/np.pi, label="yz")
    plt.legend()
    plt.ylabel("deg")
    plt.show()
