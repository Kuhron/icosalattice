# looking at paths of successive points to see how the direction changes, how far apart each pair of points is, etc.

import numpy as np
import matplotlib.pyplot as plt


def plot_distances_and_angles(distances, angles):
    plt.subplot(2, 1, 1)
    plt.plot(distances)
    plt.ylabel("distance")
    plt.subplot(2, 1, 2)
    plt.plot(angles * 180/np.pi)
    plt.ylabel("deg")
    plt.show()

