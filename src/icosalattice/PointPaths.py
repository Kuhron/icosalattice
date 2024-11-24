# functions about paths between two points
# or series of points forming a path


import numpy as np

import icosalattice.PointCodeArithmetic as pca


def get_point_path(pc_init, pc_final, direction):
    pcs = [pc_init]
    last_pc_added = pc_init
    pcs_seen = {pc_init}
    while last_pc_added != pc_final:
        next_pc = pca.add_direction_to_point_code(last_pc_added, direction)
        if next_pc in pcs_seen:
            raise RuntimeError(f"loop detected: point code {next_pc} is already in the path in the {direction} direction from {pc_init}, so the destination {pc_final} will never be reached")
        pcs.append(next_pc)
        pcs_seen.add(next_pc)
        last_pc_added = next_pc
    return pcs


def get_stepwise_path_distances_and_angles_2d(xs, ys):
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


def get_stepwise_path_distances_and_angles_3d(xs, ys, zs):
    assert len(xs) == len(ys) == len(zs)
    if len(xs) < 1:
        raise ValueError("need at least 2 points")
    dxs = np.diff(xs, n=1)
    dys = np.diff(ys, n=1)
    dzs = np.diff(zs, n=1)
    distances = (dxs**2 + dys**2 + dzs**2)**0.5
    angles_xy = np.arctan2(dys, dxs)
    angles_xz = np.arctan2(dzs, dxs)
    angles_yz = np.arctan2(dzs, dys)
    return distances, angles_xy, angles_xz, angles_yz
