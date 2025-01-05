# figuring out the correspondence between point code / peel coordinates, xyz, and location projected onto the face plane

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from functools import reduce

import icosalattice.StartingPoints as sp
import icosalattice.IcosahedronMath as icm
import icosalattice.CoordinatesByAncestry as anc
import icosalattice.PeelCoordinates as pe
import icosalattice.PointCodeArithmetic as pca
import icosalattice.MathUtil as mu
import icosalattice.GeneratePointCodes as gpc
import icosalattice.PlotPointLocations as ppl
from icosalattice.PointPaths import get_point_path, get_stepwise_path_distances_and_angles_2d
from icosalattice.PlotPaths import plot_distances_and_angles_2d
from icosalattice.CoordinatesOfPointCode import METHOD_NAME_TO_FUNCTION_POINT_CODE_TO_XYZ



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
    # plot a half-peel and some basic points
    ls = [0, 0, 0, 0.5, 0.5, 0.5, 1, 1, 1]
    ds = [0, 0.5, 1] * 3
    labels = ["C", "C3", "L", "C1", "C2", "L1", "A", "K1", "K"]
    ppl.plot_points_by_peel_coordinates_on_one_half_peel(ls, ds, labels)

    # # test if line segments are parallel on the face plane
    # reference_pair_l = ["C", "A"]
    # test_pairs_l = [
    #     ["C2", "K1"], ["C02", "K11"], ["C22", "K01"], ["C002", "K111"], ["C022", "K101"], ["C202", "K011"], ["C222", "K001"],  # both points on edges
    #     ["C13", "K11"], ["C2", "C21"],  # one point on an edge and the other inside face
    #     ["C13", "C12"], ["C012", "C112"],  # both points inside face
    # ]
    # reference_pair_dl = ["C", "K"]
    # test_pairs_dl = [
    #     ["C1", "K1"], ["C11", "K11"], ["C111", "K111"], ["C01", "K01"], ["C101", "K101"], ["C011", "K011"], ["C001", "K001"],  # both points on edges
    #     ["C1", "C12"], ["C01", "C13"], ["C12", "K1"],  # one point on an edge and the other inside face
    #     ["C21", "C212"], ["C13", "C21"],  # both points inside face
    # ]
    # reference_pair_d = ["A", "K"]  # this is only in the D direction from the perspective of the CAKL half-peel; from K's perspective it is the R direction!
    # test_pairs_d = [
    #     ["C1", "C2"], ["C01", "C02"], ["C11", "C22"], ["C001", "C002"], ["C011", "C022"], ["C101", "C202"], ["C111", "C222"],  # both points on edges
    #     ["C1", "C13"], ["C12", "C22"],  # one point on an edge and the other inside face
    #     ["C12", "C21"], ["C0103", "C0133"],  # both points inside face
    # ]
    # for reference_pair, test_pairs in zip(
    #     [reference_pair_l, reference_pair_dl, reference_pair_d],
    #     [test_pairs_l, test_pairs_dl, test_pairs_d],
    # ):
    #     pair0 = reference_pair
    #     for pair1 in test_pairs:
    #         if vectors_between_pairs_of_point_codes_are_parallel(pair0, pair1):
    #             print(f"{pair0} vs {pair1}: IS parallel")
    #         else:
    #             print(f"{pair0} vs {pair1}: is NOT parallel")

    # CHOSEN_METHOD = "ebs1"
    # CHOSEN_METHOD = "upg1"
    # CHOSEN_METHOD = "cpg1"
    CHOSEN_METHOD = "rta1"
    func_pc_to_xyz = METHOD_NAME_TO_FUNCTION_POINT_CODE_TO_XYZ[CHOSEN_METHOD]

    # observe how the points lie along lines or curves on the face
    # pcs = gpc.get_all_point_codes_from_ancestor_at_iteration(ancestor_pc="C", iterations=6)
    pcs = gpc.get_all_point_codes_on_face_at_iteration(face_name="CAKX", iterations=4, with_edges=True, with_trailing_zeros=True)
    ppl.plot_point_codes_on_half_peel_face_planes(pcs, face_name="CAKX", func_pc_to_xyz=func_pc_to_xyz, with_labels=False)

    pcs2 = reduce((lambda x,y: x+y), [gpc.get_all_point_codes_from_ancestor_at_iteration(ancestor_pc=apc, iterations=6) for apc in "CEGIK"]) + ["A"]
    ppl.plot_point_codes_on_sphere_3d(pcs2, func_pc_to_xyz=func_pc_to_xyz, with_labels=False)

    # TODO write a function that maps from an idealized face plane triangle (ACK centered on the origin)
    # - to an idealized sphere surface section (ACK where the vertices are centered on the origin and it bows upward in the z direction)
    # - and plot some points on both of those views
    # - there should only be one such mapping between the sphere surface and the face plane: the one defined by projection to/from the sphere center

    # TODO write a function that maps from the triangle to itself subject to the constraints that:
    # - it is one-to-one
    # - it is symmetric for any of the symmetries of the triangle (D3 group): 120 deg turns about the centroid, reflections
    # - things that map to themselves: the vertices, the centroid, the midpoints of the three edges
    # - the edges also map to themselves (but not necessarily pointwise, and in fact they shouldn't: there needs to be distortion between vertex and midpoint)
    # - there can be many such functions, play with them and see how the points look on the sphere surface when plotted there

    # # plot the points used in the line segment tests, so I can look at where they actually are after distortion induced by projection onto the face plane
    # pcs = sorted(set(reduce(lambda x,y: x+y, [reference_pair_l, reference_pair_dl, reference_pair_d] + test_pairs_l + test_pairs_dl + test_pairs_d)))
    # plot_point_codes_on_half_peel_face_planes(pcs, face_name="CAKX", func_pc_to_xyz=func_pc_to_xyz)
    # plot_point_codes_on_sphere_3d(pcs, with_labels=True)


    # try a "line" of points only going in one icosa direction
    pc_init, pc_final, direction = "C010100000", "K010100000", 2
    # with func_pc_to_xyz from ancestry (edge bisection), we get weird fractals such as the "alternating staircase fractal"
    # pc_init, pc_final, direction = "C010000000", "K010000000", 2  # three straight segments
    # pc_init, pc_final, direction = "C010000001", "K010000001", 2  # three straight segments with spike fractal
    # pc_init, pc_final, direction = "C02020202", "K10101011", 1  # alternating staircase fractal
    # pc_init, pc_final, direction = "C0101010101", "K0101010101", 2  # alternating staircase fractal
    # pc_init, pc_final, direction = "C0100011001", "K0100011001", 2  # mostly straight with some abrupt jumps and looking a little similar to the staircase
    pcs = get_point_path(pc_init, pc_final, direction)[:-1]
    ppl.plot_point_codes_on_half_peel_face_planes(pcs, face_name="CAKX", func_pc_to_xyz=func_pc_to_xyz, with_labels=False)
    ls, ds = pe.get_adjusted_peel_coordinates_of_point_codes_on_face(pcs, face_name="CAKX", func_pc_to_xyz=func_pc_to_xyz)
    distances, angles = get_stepwise_path_distances_and_angles_2d(ls, ds)
    plot_distances_and_angles_2d(distances, angles)
