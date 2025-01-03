import random
from functools import reduce

import icosalattice.Iterations as it
import icosalattice.PointCodeArithmetic as pca
import icosalattice.StartingPoints as sp
from icosalattice.Faces import FACE_NAMES, get_directionality_of_face, select_point_codes_on_face



def get_descendants_of_point_code_using_directions(ancestor_pc, directions, iterations, with_trailing_zeros=True):
    its_of_ancestor = it.get_iteration_number_from_point_code(ancestor_pc)
    if iterations < its_of_ancestor:
        raise ValueError(f"{ancestor_pc = } already has {its_of_ancestor} iterations, so cannot get descendants with only {iterations} iterations")
    
    iterations_needed = iterations - its_of_ancestor
    pcs = [ancestor_pc]
    for i in range(iterations_needed):
        pcs = reduce(lambda x,y: x+y, [[pc + x for x in directions] for pc in pcs])
    if not with_trailing_zeros:
        pcs = [pca.strip_trailing_zeros(pc) for pc in pcs]
    assert len(pcs) == len(set(pcs)), "shouldn't have had duplicates"
    return pcs


def get_all_point_codes_from_ancestor_at_iteration(ancestor_pc, iterations, with_trailing_zeros=True):
    if ancestor_pc in ["A", "B"]:
        if with_trailing_zeros:
            pc = ancestor_pc + "0" * iterations
        else:
            pc = ancestor_pc
        return [pc]
    return get_descendants_of_point_code_using_directions(ancestor_pc, "0123", iterations=iterations, with_trailing_zeros=with_trailing_zeros)


def get_all_point_codes_on_face_at_iteration(face_name, iterations, with_edges=True, with_trailing_zeros=True):
    directionality = get_directionality_of_face(face_name)
    p0, p1, p2, p3 = face_name

    # start off with getting the normal ancestors of the point code which are on this face
    # we are going to throw away half of them, which can be optimized out later if needed
    pcs = get_all_point_codes_from_ancestor_at_iteration(ancestor_pc=p0, iterations=iterations, with_trailing_zeros=with_trailing_zeros)
    pcs = select_point_codes_on_face(pcs, face_name)

    if with_edges:
        # add missing edges and vertices
        ring = sp.get_starting_point_ring_from_point_code(p0)
        if ring == "northern_ring":
            if directionality == "up":
                # e.g. CAKX
                # CA and CK are included in 1 and 2 directions from C
                # missing A, K, and KA (K+1)
                assert p1 == "A"
                missing_lone_vertex = p1  # A
                missing_edge_apc = p2  # K
                missing_edge_directions = "01"  # from K to A
            else:  # directionality == "down"
                # e.g. CXKL
                # CK and CL are included in 2 and 3 directions from C
                # missing K, L, and LK (L+1)
                missing_lone_vertex = p2  # K
                missing_edge_apc = p3  # L
                missing_edge_directions = "01"  # from L to K
        elif ring == "southern_ring":
            if directionality == "up":
                # e.g. DCLX
                # DC and DL are included in 1 and 2 directions from D
                # missing C, L, and CL (C+3)
                missing_lone_vertex = p2  # L
                missing_edge_apc = p1  # C
                missing_edge_directions = "03"  # from C to L
            else:  # directionality == "down"
                # e.g. DXLB
                # DL and DB are included in 2 and 3 directions from D
                # missing L, B, and LB (L+3)
                assert p3 == "B"
                missing_lone_vertex = p3  # B
                missing_edge_apc = p2  # L
                missing_edge_directions = "03"  # from L to B
        else:
            raise ValueError(f"invalid {ring = }")
        
        edge_pcs = []
        if with_trailing_zeros:
            missing_lone_vertex = pca.pad_with_trailing_zeros(missing_lone_vertex, iterations=iterations)
        edge_pcs.append(missing_lone_vertex)
        
        edge_pcs += get_descendants_of_point_code_using_directions(ancestor_pc=missing_edge_apc, directions=missing_edge_directions, iterations=iterations, with_trailing_zeros=with_trailing_zeros)
        pcs += edge_pcs
    
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
