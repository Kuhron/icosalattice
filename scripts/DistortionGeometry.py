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


def plot_peel_coordinates_on_one_half_peel(ls, ds, labels=None):
    if labels is None:
        labels = [None for x in ls]
    assert len(ls) == len(ds) == len(labels)

    plt.subplot(1, 2, 1)
    # lines for triangles
    x_mid = 0
    x_right = 1/2
    x_left = -x_right
    y_mid = 0
    y_top = np.sqrt(3)/2
    y_bottom = -y_top
    plt.plot([x_right, x_mid, x_left, x_right, x_mid, x_left], [y_mid, y_top, y_mid, y_mid, y_bottom, y_mid])
    plt.gca().set_aspect("equal")

    dxl = x_mid - x_right
    dyl = y_top - y_mid
    dxd = x_mid - x_right
    dyd = y_bottom - y_mid
    x0 = x_right
    y0 = y_mid

    # I'm sure numpy can do all this in one matrix operation but I can't be bothered right now
    xs = []
    ys = []
    for l, d, label in zip(ls, ds, labels):
        if l is None:
            assert d is None
            continue
        dx = l * dxl + d * dxd
        dy = l * dyl + d * dyd
        x = x0 + dx
        y = y0 + dy
        xs.append(x)
        ys.append(y)
        plt.gca().annotate(label, (x, y))
    plt.scatter(xs, ys)

    # view as half-peel sheared into a square
    plt.subplot(1, 2, 2)
    # lines for triangles
    x_right = 1/2
    x_left = -x_right
    y_top = 1/2
    y_bottom = -y_top
    plt.plot([x_right, x_left, x_left, x_right, x_right, x_left], [y_top, y_top, y_bottom, y_top, y_bottom, y_bottom])
    plt.gca().set_aspect("equal")

    xs = []
    ys = []
    for l, d, label in zip(ls, ds, labels):
        if l is None:
            assert d is None
            continue
        x = 1/2 - l
        y = 1/2 - d
        xs.append(x)
        ys.append(y)
        plt.gca().annotate(label, (x,y))
    plt.scatter(xs, ys)

    plt.show()


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


def get_peel_coordinates_of_point_codes_on_face(pcs, face_name):
    ls = []
    ds = []
    for pc in pcs:
        xyz = anc.get_xyz_from_point_code_using_ancestry(pc, as_array=True)
        try:
            l, d = pe.get_peel_coordinates_of_point_from_face_name(xyz, face_name, allow_one=True)
        except mu.InvalidVectorDecompositionException:
            l, d = None, None
        ls.append(l)
        ds.append(d)
    return ls, ds


def plot_point_codes_on_half_peel_face_planes(pcs, face_name, with_labels=True):
    ls, ds = get_peel_coordinates_of_point_codes_on_face(pcs, face_name=face_name)
    ls_to_plot = []
    ds_to_plot = []
    labels = []
    for pc, l, d in zip(pcs, ls, ds):
        if l is None:
            assert d is None
            continue
        ls_to_plot.append(l)
        ds_to_plot.append(d)
        if with_labels:
            labels.append(pc)
    plot_peel_coordinates_on_one_half_peel(ls, ds, labels=labels if with_labels else None)


def plot_point_codes_on_sphere_3d(pcs, with_labels=True):
    xyzs = [anc.get_xyz_from_point_code_using_ancestry(pc, as_array=False) for pc in pcs]
    xs, ys, zs = zip(*xyzs)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.scatter(xs, ys, zs)
    ax.set(
        xticklabels=[],
        yticklabels=[],
        zticklabels=[],
    )
    if with_labels:
        for pc, x, y, z in zip(pcs, xs, ys, zs):
            ax.text(x, y, z, pc)
    plt.show()


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
    plot_peel_coordinates_on_one_half_peel(ls, ds, labels)

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
    pcs = ["C"]
    for i in range(6):
        pcs = reduce(lambda x,y: x+y, [[pc + x for x in "0123"] for pc in pcs])
    pcs = sorted(set(pca.strip_trailing_zeros(pc) for pc in pcs))
    plot_point_codes_on_half_peel_face_planes(pcs, face_name="CAKX", with_labels=False)
    plot_point_codes_on_sphere_3d(pcs, with_labels=False)

    # try a "line" of points only going in one icosa direction
    pc_init, pc_final, direction = "C101010101", "C202020202", 3  # alternating staircase fractal
    # pc_init, pc_final, direction = "C010000000", "K010000000", 2  # three straight segments
    # pc_init, pc_final, direction = "C010000001", "K010000001", 2  # three straight segments with spike fractal
    # pc_init, pc_final, direction = "C02020202", "K10101011", 1  # alternating staircase fractal
    # pc_init, pc_final, direction = "C0101010101", "K0101010101", 2  # alternating staircase fractal
    # pc_init, pc_final, direction = "C0100011001", "K0100011001", 2  # mostly straight with some abrupt jumps and looking a little similar to the staircase
    pcs = [pc_init]
    while pcs[-1] != pc_final:
        pcs.append(pca.add_direction_to_point_code(pcs[-1], direction))
    plot_point_codes_on_half_peel_face_planes(pcs, face_name="CAKX", with_labels=False)
    ls, ds = get_peel_coordinates_of_point_codes_on_face(pcs, face_name="CAKX")
    distances, angles = get_stepwise_path_distances_and_angles(ls, ds)
    plt.subplot(2, 1, 1)
    plt.plot(distances)
    plt.ylabel("distance")
    plt.subplot(2, 1, 2)
    plt.plot(angles * 180/np.pi)
    plt.ylabel("deg")
    plt.show()


    # name of original method for locating child points: Bisection From Parents (or maybe just "Bisection" for short)

    # TODO one idea for redoing where the child points are:
    # - method name: Corrected Plane Gridding
    # - use peel coordinates as primary
    # - but dilate where the l and d proportions are by using the trig formula for lp (so the same delta l moves less near the center and more near the corners)
    # - draw triangular grid on the face plane based solely on these lp points on the three edges, so the grid lines are straight on the face plane
    # - project the resulting points to the sphere
    # - see what they look like and how far apart they are, hopefully the lp adjustment makes them relatively uniformly spaced at least within a face
    # - and hopefully the lines put the points in neat rows on the sphere surface, and most or all refraction occurs at actual face boundaries
    # - (want no Sierpinski artifacts in the distribution of points within a face)

    # TODO another idea for child point locations:
    # - method name: Arc Gridding
    # - divide the sphere edges into n equal segments based on which iteration you are at (e.g. 8 segments)
    # - then connect corresponding points along these edges using great circle paths
    # - hopefully the gc paths meet at sixfold vertices
    # - and I suspect this method is equivalent to Grid Distortion

    # the method of doing the gridding on the face plane uniformly and then projecting that out onto the sphere, introducing distortion,
    # - could be called Uncorrected Plane Gridding

    