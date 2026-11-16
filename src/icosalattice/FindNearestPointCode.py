import numpy as np

import icosalattice.MapCoordinateMath as mcm
from icosalattice.Iterations import get_iteration_number_from_point_code, adjust_point_code_to_iteration
from icosalattice.Adjacency import get_neighbors_of_point_code
from icosalattice.CoordinatesOfPointCode import get_xyz_from_point_code, get_xyzs_from_point_codes, METHOD_NAMES
from icosalattice.StartingPoints import STARTING_POINT_CODES
from icosalattice.RandomUtil import get_pseudorandom_point_codes



def get_point_code_from_xyz(xyz, pc_to_xyz_method):
    raise NotImplementedError("do some kind of binary search thing")


def get_point_code_from_latlon(latlon, pc_to_xyz_method):
    xyz = mcm.get_xyz_from_latlon(latlon)
    return get_point_code_from_xyz(xyz, pc_to_xyz_method)


# # TODO test/refactor; this is from old IcosahedronMath.py
# def get_nearest_icosa_point_at_lower_iteration(pc, iterations, xyzg):
#     pc_iter = get_iteration_number_from_point_code(pc)
#     if iterations == pc_iter:
#         return pc
#     elif iterations > pc_iter:
#         return pc + "0"*(iterations - pc_iter)
#     else:
#         ancestor = pc[:iterations+1]
#         adj = get_adjacency_from_point_code(ancestor)
#         pcs = adj + [ancestor]
#         # find the closest of these, it shouldn't be anything else (I think, haven't proven that)
#         xyzs = get_xyzs_from_point_codes(pcs, xyzg)
#         xyz = get_xyz_from_point_code(pc, xyzg)
#         ds = [mcm.xyz_distance(xyz, xyz1) for xyz1 in xyzs]
#         index = np.argmin(ds)
#         return pcs[index]


def get_nearest_point_code_to_latlon(latlon, maximum_distance, planet_radius, pc_to_xyz_method):
    xyz = mcm.get_xyz_from_latlon(latlon)
    return get_nearest_point_code_to_xyz(xyz, maximum_distance, planet_radius, pc_to_xyz_method)


# TODO test/refactor; this is from old IcosahedronMath.py
def get_nearest_point_code_to_xyz(xyz, maximum_distance, planet_radius, pc_to_xyz_method):
    # print("getting nearest icosa point to {}".format(xyz))
    max_distance_normalized = maximum_distance / planet_radius
    candidate_pcs = STARTING_POINT_CODES
    iteration = 0
    while True:
        candidate_pcs = [adjust_point_code_to_iteration(pc, iteration) for pc in candidate_pcs]
        # print(f"{iteration = }; {candidate_pcs = }")
        nearest_candidate_pc, distance_normalized = get_nearest_neighbor_to_xyz(xyz, candidate_pcs, pc_to_xyz_method)
        # print("nearest candidate is {} at distance of {}".format(nearest_candidate_pc, distance_normalized))
        if distance_normalized <= max_distance_normalized:
            # print("-- done getting nearest icosa point to {}".format(xyz))
            distance_in_units = distance_normalized * planet_radius
            return nearest_candidate_pc, distance_normalized, distance_in_units

        iteration += 1
        if iteration > 30:
            raise RuntimeError("while loop ran too many times")
        nearest_candidate_neighbor_pcs = get_neighbors_of_point_code(nearest_candidate_pc)
        nearest_candidate_neighbor_pcs = [x for x in nearest_candidate_neighbor_pcs if x is not None]
        candidate_pcs = nearest_candidate_neighbor_pcs + [nearest_candidate_pc]


# TODO test/refactor; this is from old IcosahedronMath.py
def get_nearest_neighbor_to_latlon(latlon, candidate_pcs, pc_to_xyz_method):
    xyz = mcm.get_xyz_from_latlon(latlon)
    return get_nearest_neighbor_to_xyz(xyz, candidate_pcs, pc_to_xyz_method)


# TODO test/refactor; this is from old IcosahedronMath.py
def get_nearest_neighbor_to_xyz(xyz, candidate_pcs, pc_to_xyz_method):
    candidates_dict = {tuple(get_xyz_from_point_code(pc, pc_to_xyz_method, as_array=False)): pc for pc in candidate_pcs}
    candidates_xyz = list(candidates_dict.keys())
    nn_xyz, d = get_nearest_neighbor_xyz_to_xyz(xyz, candidates_xyz)
    return candidates_dict[nn_xyz], d


# TODO test/refactor; this is from old IcosahedronMath.py
def get_nearest_neighbor_xyz_to_xyz(xyz, candidates_xyz):
    p = np.array(xyz)
    assert p.shape == (3,), p.shape
    ps = np.array(candidates_xyz)
    n = len(candidates_xyz)
    if n == 0:
        raise ValueError("empty candidates list")
    assert ps.shape == (n, 3), f"ps should have shape ({n}, 3) but got {ps.shape}:\n{ps}"
    dx = p-ps
    d = sum((dx**2).T) ** 0.5
    assert d.shape == (n,)
    nn_index = np.argmin(d)
    min_d = min(d)
    nn_arr = ps[nn_index]
    nn_xyz = tuple(nn_arr)
    return nn_xyz, min_d



if __name__ == "__main__":
    latlons = [
        [47.37444, 8.54111],  # Zurich, Switzerland
        [10.77556, 106.70194],  # Ho Chi Minh City, Vietnam
        [-54.807222, -68.308056],  # Ushuaia, Argentina
        [44.052222, -123.0925],  # Eugene, Oregon, United States
        [-18.141600, 178.441895],  # Suva, Fiji
        [90, 49.3338],  # North Pole with irrelevant longitude
        [-90, -71.64665],  # South Pole with irrelevant longitude
        [-89.9208, -71.64665],  # close to South Pole
    ]

    for latlon in latlons:
        print(f"{latlon = }")
        for method in METHOD_NAMES:
            # pc = get_point_code_from_latlon(latlon, method=method)
            nearest_candidate_pc, distance_normalized, distance_in_units = get_nearest_point_code_to_latlon(latlon, maximum_distance=1, planet_radius=6371, pc_to_xyz_method=method)
            print(f"  {method} chose pc {nearest_candidate_pc} at distance {distance_in_units:.3f} km ({distance_normalized:.6f} on unit sphere)")

    # FIXME bugs:
    # - all found point codes end in 0
    