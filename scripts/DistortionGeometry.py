# figuring out the correspondence between point code / peel coordinates, xyz, and location projected onto the face plane
# compare different methods of placing points onto face planes based on their raw peel coordinates (trying to get nice distribution of lattice on the sphere)

import numpy as np

from icosalattice.CoordinatesOfPointCode import METHOD_NAME_TO_FUNCTION_POINT_CODE_TO_XYZ
import icosalattice.EvaluatePointPlacementMethods as epp
import icosalattice.Adjacency as adj



def get_vector_between_xyzs_from_point_codes(pc0, pc1, func_pc_to_xyz):
    xyz0 = func_pc_to_xyz(pc0, as_array=True)
    xyz1 = func_pc_to_xyz(pc1, as_array=True)
    return xyz1 - xyz0


def vectors_between_pairs_of_point_codes_are_parallel(pair0, pair1, func_pc_to_xyz):
    pc00, pc01 = pair0
    pc10, pc11 = pair1
    v0 = get_vector_between_xyzs_from_point_codes(pc00, pc01, func_pc_to_xyz)
    v1 = get_vector_between_xyzs_from_point_codes(pc10, pc11, func_pc_to_xyz)
    cross = np.cross(v0, v1)
    mag = np.linalg.norm(cross)
    return mag < 1e-9



if __name__ == "__main__":
    # # plot a half-peel and some basic points
    # ls = [0, 0, 0, 0.5, 0.5, 0.5, 1, 1, 1]
    # ds = [0, 0.5, 1] * 3
    # labels = ["C", "C3", "L", "C1", "C2", "L1", "A", "K1", "K"]
    # ppl.plot_points_by_peel_coordinates_on_one_half_peel(ls, ds, labels)

    # CHOSEN_METHOD = "ebs1"
    # CHOSEN_METHOD = "upg1"
    CHOSEN_METHOD = "cpg1"
    # CHOSEN_METHOD = "rta1"
    func_pc_to_xyz = METHOD_NAME_TO_FUNCTION_POINT_CODE_TO_XYZ[CHOSEN_METHOD]
    
    # epp.plot_point_distribution_on_one_face(func_pc_to_xyz, iterations=6)

    # epp.plot_point_distribution_on_northern_half_peels(func_pc_to_xyz, iterations=6)

    # try a "line" of points only going in one icosa direction
    pc_init, pc_final, direction = "C010100000", "K010100000", 2
    # with func_pc_to_xyz from ancestry (edge bisection), we get weird fractals such as the "alternating staircase fractal"
    # pc_init, pc_final, direction = "C010000000", "K010000000", 2  # three straight segments
    # pc_init, pc_final, direction = "C010000001", "K010000001", 2  # three straight segments with spike fractal
    # pc_init, pc_final, direction = "C02020202", "K10101011", 1  # alternating staircase fractal
    # pc_init, pc_final, direction = "C0101010101", "K0101010101", 2  # alternating staircase fractal
    # pc_init, pc_final, direction = "C0100011001", "K0100011001", 2  # mostly straight with some abrupt jumps and looking a little similar to the staircase
    # epp.plot_trajectory_between_two_points_on_face_plane(pc_init, pc_final, direction, func_pc_to_xyz)

    epp.report_neighbor_angle_and_distance_statistics(func_pc_to_xyz, iterations=6)
