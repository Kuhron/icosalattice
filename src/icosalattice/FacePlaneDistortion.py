# functions about the distortion introduced when projecting from the plane containing an icosahedron face (whose vertices are on the sphere surface)
# or when projecting from this "face plane" back onto the sphere surface


import numpy as np
import matplotlib.pyplot as plt

import icosalattice.IcosahedronMath as icm
import icosalattice.StartingPoints as sp


ALPHA = icm.ANGLE_BETWEEN_VERTICES_RAD
H = np.cos(ALPHA/2)  # height from center of sphere to middle of edge of face plane
W = 2*np.sin(ALPHA/2)  # length of edge of face plane


def get_x_from_theta(theta):
    # `x`: distance from the middle of the face plane edge to the projected point
    return H * np.tan(ALPHA/2 - theta)


def get_theta_from_x(x):
    return ALPHA/2 - np.atan2(x/H)


def get_lp_from_theta(theta):
    # `lp`: length along face plane edge (from one vertex) of the projection of a point above it on the sphere
    return W/2 - get_x_from_theta(theta)


def get_theta_from_lp(lp):
    # `theta`: angle (from perspective of sphere's center) from vertex to point at lp along the face plane's edge
    x = W/2 - lp
    return get_theta_from_x(x)


def get_lp_from_theta_proportion(a):
    return get_lp_from_theta(a * ALPHA)


def get_theta_proportion_from_lp(lp):
    theta = get_theta_from_lp(lp)
    return theta / ALPHA


def get_lp_proportion_from_theta_proportion(a):
    return get_lp_from_theta(a * ALPHA) / W


def get_theta_proportion_from_lp_proportion(a):
    return get_theta_proportion_from_lp(a * W)



if __name__ == "__main__":
    p_a = sp.STARTING_POINTS[0]
    p_c = sp.STARTING_POINTS[2]
    p_k = sp.STARTING_POINTS[10]
    xyz_a = p_a.xyz(as_array=True)
    xyz_c = p_c.xyz(as_array=True)
    xyz_k = p_k.xyz(as_array=True)

    a_dot_c = np.dot(xyz_a, xyz_c)
    alpha_by_dot_product = np.acos(a_dot_c) / (1 * 1)
    assert ALPHA == alpha_by_dot_product

    h_from_radicals = np.sqrt(1/10 * (5 + np.sqrt(5)))
    assert H == h_from_radicals

    assert np.isclose(W, np.linalg.norm(xyz_a - xyz_c), atol=1e-9)
    assert np.isclose(W, np.linalg.norm(xyz_a - xyz_k), atol=1e-9)
    assert np.isclose(W, np.linalg.norm(xyz_c - xyz_k), atol=1e-9)

    assert H**2 + (W/2)**2 == 1**2 # Pythagoras for triangle (sphere center, vertex, middle of face plane edge)

    # plot how lp deviates from assuming linear movement along face plane edge
    proportions = np.linspace(0, 1, 101)
    lp_proportions = get_lp_proportion_from_theta_proportion(proportions)
    deviations_from_linear = lp_proportions - proportions
    assert deviations_from_linear[0] == deviations_from_linear[50] == deviations_from_linear[100] == 0
    plt.plot(proportions, deviations_from_linear)
    plt.show()
