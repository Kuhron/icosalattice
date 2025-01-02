# functions about the distortion introduced when projecting from the plane containing an icosahedron face (whose vertices are on the sphere surface)
# or when projecting from this "face plane" back onto the sphere surface


import numpy as np
import matplotlib.pyplot as plt

import icosalattice.IcosahedronMath as icm
import icosalattice.StartingPoints as sp


ALPHA = icm.ANGLE_BETWEEN_VERTICES_RAD

# H, the height from the center of the sphere to the middle of an edge of a face plane
H = np.cos(ALPHA/2)
assert np.isclose(H, (1/10 * (5 + 5**0.5))**0.5, atol=1e-9)  # found from putting decimal in Wolfram Alpha

# W, the length of the edge of the face plane
# = length of the straight line in 3D between two neighboring icosa vertices
W = 2*np.sin(ALPHA/2)  # length of edge of face plane
assert np.isclose(W, (2/5 * (5 - 5**0.5))**0.5, atol=1e-9)  # found from putting decimal in Wolfram Alpha

# B, the length from a vertex to the middle of the opposite edge of a face plane
B = (1/2 * 3**0.5) * W
assert np.isclose(B, (3/10 * (5 - 5**0.5))**0.5, atol=1e-9)  # found from putting decimal in Wolfram Alpha
# for WA: b = sqrt(3/10 * (5 - sqrt(5)))

# G, the length from the center of the sphere to the centroid of a face plane
# G = ((-(B**4) + 2*(B**2)*(H**2+1) - (H**2-1)**2) ** 0.5) / (2*B)

G = (1 - ((B**2 - H**2 + 1)**2)/(4 * B**2))**0.5
assert np.isclose(G, (1/15 * (5 + 2*5**0.5))**0.5, atol=1e-9)  # found from iteratively simplifying this expression on paper and in Wolfram Alpha

# for WA: g = sqrt(1/15 * (5 + 2*sqrt(5)))
#         g^2 = 1/15 * (5 + 2*sqrt(5))

# check correctness of triangles (sphere center, A, face plane centroid) and (sphere center, mid(CK), face plane centroid)
b1_from_wolfram = (B**2 - H**2 + 1)/(2*B)
b1_from_triangle = (1 - G**2)**0.5
b2 = B - b1_from_triangle
assert np.isclose(b1_from_wolfram, b1_from_triangle, atol=1e-9)
assert np.isclose(G**2 + b1_from_wolfram**2, 1, atol=1e-9)
assert np.isclose(G**2 + b2**2, H**2, atol=1e-9)
del b1_from_wolfram
del b1_from_triangle
del b2

assert np.isclose(B/G, 1/2*(9-3*5**0.5), atol=1e-9)
# for WA: b = g/2*(9-3*sqrt(5))

# GAMMA, the angle from a corner of the face plane to its center (angle from perspective of the sphere center)
GAMMA = np.arccos(G/1)
assert np.isclose(GAMMA, np.arctan(3 - 5**0.5), atol=1e-9)  # found from putting decimal in Wolfram Alpha
# for WA: gamma = arctan(3-sqrt(5))


def get_x_from_theta(theta):
    # `x`: distance from the middle of the face plane edge to the projected point
    return H * np.tan(ALPHA/2 - theta)


def get_theta_from_x(x):
    return ALPHA/2 - np.atan2(x, H)


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
    return get_lp_from_theta_proportion(a) / W


def get_theta_proportion_from_lp_proportion(a):
    return get_theta_proportion_from_lp(a * W)



if __name__ == "__main__":
    print(f"{ALPHA = }")
    print(f"{W = }")
    print(f"{H = }")
    print(f"{B = }")
    print(f"{G = }")
    print(f"{GAMMA = }")


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
