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
from icosalattice.PointPaths import get_point_path
from icosalattice.PlotPaths import plot_distances_and_angles



def get_vector_between_xyzs_from_point_codes(pc0, pc1):
    xyz0 = anc.get_xyz_from_point_code_using_ancestry(pc0, as_array=True)
    xyz1 = anc.get_xyz_from_point_code_using_ancestry(pc1, as_array=True)
    return xyz1 - xyz0


def vectors_between_pairs_of_point_codes_are_parallel(pair0, pair1):
    pc00, pc01 = pair0
    pc10, pc11 = pair1
    v0 = get_vector_between_xyzs_from_point_codes(pc00, pc01)
    v1 = get_vector_between_xyzs_from_point_codes(pc10, pc11)
    cross = np.cross(v0, v1)
    mag = np.linalg.norm(cross)
    return mag < 1e-9


def get_stepwise_path_distances_and_angles(xs, ys):
    # see how a path of points looks
    # for qualitatively investigating what is going on with what should be lines or curves of points along a single direction path on a face plane
    # (e.g., D from C100 to C200)
    assert len(xs) == len(ys)
    if len(xs) < 1:
        raise ValueError("need at least 2 points")
    dxs = np.diff(xs, n=1)
    dys = np.diff(ys, n=1)
    distances = (dxs**2 + dys**2)**0.5
    angles = np.arctan2(dys, dxs)
    return distances, angles



if __name__ == "__main__":
    # angle between two vertices of icosahedron
    alpha = icm.ANGLE_BETWEEN_VERTICES_RAD

    p_a = sp.STARTING_POINTS[0]
    p_c = sp.STARTING_POINTS[2]
    p_k = sp.STARTING_POINTS[10]
    xyz_a = p_a.xyz(as_array=True)
    xyz_c = p_c.xyz(as_array=True)
    xyz_k = p_k.xyz(as_array=True)

    a_dot_c = np.dot(xyz_a, xyz_c)
    alpha_by_dot_product = np.acos(a_dot_c) / (1 * 1)
    assert alpha == alpha_by_dot_product

    h = np.cos(alpha/2)  # height from center of sphere to middle of edge of face plane
    h_from_radicals = np.sqrt(1/10 * (5 + np.sqrt(5)))
    assert h == h_from_radicals

    w = 2*np.sin(alpha/2)  # length of edge of face plane
    assert np.isclose(w, np.linalg.norm(xyz_a - xyz_c), atol=1e-9)
    assert np.isclose(w, np.linalg.norm(xyz_a - xyz_k), atol=1e-9)
    assert np.isclose(w, np.linalg.norm(xyz_c - xyz_k), atol=1e-9)

    assert h**2 + (w/2)**2 == 1**2 # Pythagoras for triangle (sphere center, vertex, middle of face plane edge)

    # lp is the length along face plane edge (from one vertex) of the projection of a point above it on the sphere
    get_lp = lambda theta: w/2 - get_x(theta)

    # x is the distance from the middle of the face plane edge to the projected point
    get_x = lambda theta: h * np.tan(alpha/2 - theta)

    # # plot how lp deviates from assuming linear movement along face plane edge
    # thetas = np.linspace(0, alpha, 101)
    # lps = get_lp(thetas)
    # proportions = thetas/alpha
    # deviations_from_linear = lps - (proportions * w)
    # assert deviations_from_linear[0] == deviations_from_linear[50] == deviations_from_linear[100] == 0
    # plt.plot(thetas, deviations_from_linear)
    # plt.show()

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
    
    # # plot the points used in the line segment tests, so I can look at where they actually are after distortion induced by projection onto the face plane
    # pcs = sorted(set(reduce(lambda x,y: x+y, [reference_pair_l, reference_pair_dl, reference_pair_d] + test_pairs_l + test_pairs_dl + test_pairs_d)))
    # plot_point_codes_on_half_peel_face_planes(pcs, face_name="CAKX")
    # plot_point_codes_on_sphere_3d(pcs, with_labels=True)

    # observe how the points lie along lines or curves on the face
    pcs = gpc.get_all_point_codes_from_ancestor_at_iteration(ancestor_pc="C", iterations=6)
    ppl.plot_point_codes_on_half_peel_face_planes(pcs, face_name="CAKX", with_labels=False)
    ppl.plot_point_codes_on_sphere_3d(pcs, with_labels=False)

    # try a "line" of points only going in one icosa direction
    pc_init, pc_final, direction = "C101010101", "C202020202", 3  # alternating staircase fractal
    # pc_init, pc_final, direction = "C010000000", "K010000000", 2  # three straight segments
    # pc_init, pc_final, direction = "C010000001", "K010000001", 2  # three straight segments with spike fractal
    # pc_init, pc_final, direction = "C02020202", "K10101011", 1  # alternating staircase fractal
    # pc_init, pc_final, direction = "C0101010101", "K0101010101", 2  # alternating staircase fractal
    # pc_init, pc_final, direction = "C0100011001", "K0100011001", 2  # mostly straight with some abrupt jumps and looking a little similar to the staircase
    pcs = get_point_path(pc_init, pc_final, direction)
    ppl.plot_point_codes_on_half_peel_face_planes(pcs, face_name="CAKX", with_labels=False)
    ls, ds = pe.get_peel_coordinates_of_point_codes_on_face(pcs, face_name="CAKX")
    distances, angles = get_stepwise_path_distances_and_angles(ls, ds)
    plot_distances_and_angles(distances, angles)
