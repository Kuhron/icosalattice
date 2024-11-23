import numpy as np
import functools

import icosalattice.MapCoordinateMath as mcm
import icosalattice.StartingPoints as sp
from icosalattice.Adjacency import get_adjacency_from_point_code


@functools.lru_cache(maxsize=100000)
def get_xyz_from_point_code_using_ancestry(pc, as_array=True):
    if len(pc) == 1:
        res = get_xyz_of_initial_point_code(pc)
    else:
        res = get_xyz_from_point_code_recursive(pc)

    # I don't feel like adding an "as_array" kwarg to so many functions this calls, so just going to do it here right now
    if as_array:
        return np.array(res)
    else:
        return res


def get_xyz_of_initial_point_code(pc):
    p = sp.STARTING_POINTS[sp.STARTING_POINT_CODES.index(pc)]
    return p.xyz()


def get_xyz_from_point_code_recursive(pc):
    assert len(pc) > 1, "should get latlon/xyz of initial point directly"
    res = get_xyz_of_point_code_using_parents(pc)
    return res


def get_xyz_of_point_code_using_parents(pc):
    assert len(pc) > 1, "should get latlon/xyz of initial point directly"
    xyz0, xyz1 = get_parent_xyzs_from_point_code(pc)
    return get_xyz_of_child_from_parent_xyzs(xyz0, xyz1)


def get_parent_xyzs_from_point_code(pc):
    # print(f"getting parent xyzs for {pc=}")
    p0, p1 = get_parents_from_point_code(pc)
    xyz0 = get_xyz_from_point_code_using_ancestry(p0)
    xyz1 = get_xyz_from_point_code_using_ancestry(p1)
    # print(f"\n{pc=} has parents {p0} with {xyz0=} and {p1} with {xyz1=}")
    return xyz0, xyz1


def get_xyz_of_child_from_parent_xyzs(xyz0, xyz1):
    # reduce use of UnitSpherePoint objects where they are unnecessary
    # also reduce use of latlon and conversion to/from it where it is unnecessary
    res = mcm.get_unit_sphere_midpoint_from_xyz(xyz0, xyz1)
    # print(f"midpoint of\n{xyz0=} and\n{xyz1} is\n{res}")
    return res


def get_child_index_from_point_code(pc):
    # index that the point is for its parent in the generation it was born
    # if it was the Left child, 0; if it was the DownLeft child, 1; if it was the Down child, 2
    return int(pc[-1])


def get_parents_from_point_code(pc):
    if len(pc) == 1:
        par = None
        dpar = None
    elif pc[-1] == "0":
        par = pc[:-1]
        dpar = pc[:-1]
    else:
        adj = get_adjacency_from_point_code(pc)
        ci = get_child_index_from_point_code(pc)
        par = adj[-ci]
        dpar = adj[ci]
        assert par[-1] == "0", par
        assert dpar[-1] == "0", dpar
        par = par[:-1]
        dpar = dpar[:-1]
    # print(f"parents of {pc} are {par=}, {dpar=}")
    return [par, dpar]

    # old way
    # p0 = get_parent_from_point_code(pc)
    # p1 = get_directional_parent_from_point_code(pc)
    # return [p0, p1]