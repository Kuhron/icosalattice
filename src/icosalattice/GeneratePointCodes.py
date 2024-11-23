import random
from functools import reduce

import icosalattice.Iterations as it
import icosalattice.PointCodeArithmetic as pca
import icosalattice.StartingPoints as sp



def get_all_point_codes_from_ancestor_at_iteration(ancestor_pc, iterations, with_trailing_zeros=True):
    if ancestor_pc in ["A", "B"]:
        if with_trailing_zeros:
            pc = ancestor_pc + "0" * iterations
        else:
            pc = ancestor_pc
        return [pc]
    
    its_of_ancestor = it.get_iteration_number_from_point_code(ancestor_pc)
    if iterations < its_of_ancestor:
        raise ValueError(f"{ancestor_pc = } already has {its_of_ancestor} iterations, so cannot get descendants with only {iterations} iterations")
    iterations_needed = iterations - its_of_ancestor
    pcs = [ancestor_pc]
    for i in range(iterations_needed):
        pcs = reduce(lambda x,y: x+y, [[pc + x for x in "0123"] for pc in pcs])
    if not with_trailing_zeros:
        pcs = [pca.strip_trailing_zeros(pc) for pc in pcs]
    assert len(pcs) == len(set(pcs)), "shouldn't have had duplicates"
    return pcs


def get_all_point_codes_at_iteration(iterations, with_trailing_zeros=True):
    pcs = []
    for spc in sp.STARTING_POINT_CODES:
        pcs += get_all_point_codes_from_ancestor_at_iteration(ancestor_pc=spc, iterations=iterations, with_trailing_zeros=with_trailing_zeros)
    return pcs


def get_random_point_code(min_iterations, expected_iterations, max_iterations, prefix=""):
    assert min_iterations <= expected_iterations <= max_iterations
    
    if len(prefix) == 0:
        # start the string off with the minimum iterations needed
        s = random.choice("CDEFGHIJKL")
    else:
        assert prefix[0] in "CDEFGHIJKL"
        s = prefix

    for i in range(min_iterations - len(s)):  # if negative, loop won't run (won't error)
        s += random.choice("0123")
    backup_digit = random.choice("123")  # prevent trailing zeros from taking s back below min_iterations

    iterations_used = it.get_iteration_number_from_point_code(s)
    if iterations_used < expected_iterations:
        while True:
            s += random.choice("0123")
            iterations_used = it.get_iteration_number_from_point_code(s)
            if iterations_used >= max_iterations:
                break
            if random.random() < 1/(expected_iterations - min_iterations):
                break

    # based on the number of points in this iteration, maybe replace with one of the poles
    iterations_used = it.get_iteration_number_from_point_code(s)
    n_points_this_iteration = it.get_n_points_from_iterations(iterations_used)
    prob_of_pole = 2 / n_points_this_iteration
    if min_iterations == 0 and random.random() < prob_of_pole:
        return random.choice(["A", "B"])
    
    iteration_born = it.get_iteration_born_from_point_code(pca.strip_trailing_zeros(s))
    if iteration_born < min_iterations:
        # replace the last zero with backup_digit
        s = s[:-1] + backup_digit
    else:
        s = pca.strip_trailing_zeros(s)

    return s


if __name__ == "__main__":
    for i in range(0, 5):
        pcs = get_all_point_codes_at_iteration(iterations=i, with_trailing_zeros=False)
        print(len(pcs))
    pcs = get_all_point_codes_at_iteration(iterations=2, with_trailing_zeros=True)
    print(pcs)
