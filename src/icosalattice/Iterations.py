# functions about which iteration of the lattice's generation a point is from
# note that trailing zeros on the point code can be used to indicate that it has "not moved yet" (rather than going in the 1/2/3 directions to become a child)


import numpy as np

import icosalattice.PointCodeArithmetic as pca
import icosalattice.IcosahedronMath as icm



def get_iteration_born_from_point_code(pc):
    pca.enforce_no_trailing_zeros(pc)
    return get_iteration_number_from_point_code(pc)


def get_iteration_number_from_point_code(pc):
    # this INCLUDES trailing zeros, we want to know what iteration something is at
    # for iteration born, use get_iteration_born instead
    # the initial points have length 1, just their letter A-L
    return len(pc) - 1


def get_n_points_from_iterations(n_iters: int) -> int:
    return get_exact_n_points_from_iterations(n_iters)


def get_exact_n_points_from_iterations(n_iters: int) -> int:
    return 2 + 10 * (4 ** n_iters)


def get_exact_iterations_from_n_points(n_points):
    return np.log((n_points - 2)/10) / np.log(4)  # just use change of base since np.log(arr, b) doesn't like arrays
