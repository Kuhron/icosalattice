import numpy as np
import matplotlib.pyplot as plt

import icosalattice.PeelCoordinates as pe
import icosalattice.CoordinatesByAncestry as anc



def plot_point_codes_on_half_peel_face_planes(pcs, face_name, func_pc_to_xyz, with_labels=True):
    ls, ds = pe.get_peel_coordinates_of_point_codes_on_face(pcs, face_name=face_name, func_pc_to_xyz=func_pc_to_xyz)
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
    plot_points_by_peel_coordinates_on_one_half_peel(ls, ds, labels=labels if with_labels else None)


def plot_point_codes_on_sphere_3d(pcs, with_labels=True):
    xyzs = [anc.get_xyz_from_point_code_using_ancestry(pc, as_array=False) for pc in pcs]
    plot_xyzs_on_sphere_3d(xyzs, labels=pcs if with_labels else None)


def plot_xyzs_on_sphere_3d(xyzs, labels=None):
    xs, ys, zs = zip(*xyzs)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    # ax.set_aspect("equal")  # works for 2d but not 3d
    ax.set_box_aspect([1.0, 1.0, 1.0])  # this sets aspects equal for 3d
    ax.set_xlim3d(-1.1, 1.1)
    ax.set_ylim3d(-1.1, 1.1)
    ax.set_zlim3d(-1.1, 1.1)
    ax.scatter(xs, ys, zs)
    ax.set(
        xticklabels=[],
        yticklabels=[],
        zticklabels=[],
    )
    if labels is not None:
        for label, x, y, z in zip(labels, xs, ys, zs):
            ax.text(x, y, z, label)
    plt.show()


def plot_points_by_peel_coordinates_on_one_half_peel(ls, ds, labels=None):
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


