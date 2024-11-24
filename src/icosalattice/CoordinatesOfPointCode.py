# different methods for converting a point code into coordinates
# Warning: these different methods may generate different lattices!

# name of original method for locating child points: Edge Bisection

# one idea for redoing where the child points are:
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
# - and I suspect/hope this method is equivalent to Corrected Plane Gridding


import numpy as np
import matplotlib.pyplot as plt

from icosalattice.CoordinatesByAncestry import get_xyz_from_point_code_using_ancestry
from icosalattice.GeneratePointCodes import get_all_point_codes_from_ancestor_at_iteration, get_all_point_codes_at_iteration
import icosalattice.PeelCoordinates as pe
import icosalattice.FacePlaneDistortion as distort
import icosalattice.PlotPointLocations as ppl
from icosalattice.Adjacency import get_adjacency_from_point_code
import icosalattice.MapCoordinateMath as mcm
from icosalattice.PlotPaths import plot_distances_and_angles_2d, plot_distances_and_angles_3d
from icosalattice.PointPaths import get_point_path, get_stepwise_path_distances_and_angles_2d, get_stepwise_path_distances_and_angles_3d


def get_xyz_from_point_code_using_corrected_plane_gridding(pc, as_array=True):
    # modify the l and d coordinates to try to offset the distortion introduced when projecting from face plane back onto sphere surface
    spc, l_raw, d_raw = pe.get_peel_coordinates_from_point_code(pc)
    if l_raw > d_raw:
        # on upward-pointing face
        l_modified, d_modified = adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_upward_face(l_raw, d_raw)
    elif d_raw > l_raw:
        # on downward-pointing face
        l_modified, d_modified = adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_downward_face(l_raw, d_raw)
    else:
        # on direction-2 edge, check both faces
        l_modified_0, d_modified_0 = adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_upward_face(l_raw, d_raw)
        l_modified_1, d_modified_1 = adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_downward_face(l_raw, d_raw)
        assert abs(l_modified_0 - l_modified_1) < 1e-9, f"{l_modified_0} != {l_modified_1}"
        assert abs(d_modified_0 - d_modified_1) < 1e-9, f"{d_modified_0} != {d_modified_1}"
        l_modified = (l_modified_0 + l_modified_1) / 2
        d_modified = (d_modified_0 + d_modified_1) / 2

    return pe.get_xyz_from_peel_coordinates(spc, l_modified, d_modified, as_array=as_array)


def adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_upward_face(l_raw, d_raw):
    # D vector is not in the face plane! so need in terms of L and DL

    # convert (l, d) to (l, dl) for distorting
    new_dl_raw = d_raw
    new_l_raw = l_raw - d_raw
    dl_modified = distort.get_lp_proportion_from_theta_proportion(new_dl_raw)
    l_modified = distort.get_lp_proportion_from_theta_proportion(new_l_raw)

    # now convert back to (l, d) for passing to function in PeelCoordinates module
    d_modified = dl_modified
    l_modified = l_modified + dl_modified

    return l_modified, d_modified


def adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_downward_face(l_raw, d_raw):
    # L vector is not in the face plane! so need in terms of DL and D
    
    # convert (l, d) to (dl, d) for distorting
    new_dl_raw = l_raw
    new_d_raw = d_raw - l_raw
    dl_modified = distort.get_lp_proportion_from_theta_proportion(new_dl_raw)
    d_modified = distort.get_lp_proportion_from_theta_proportion(new_d_raw)

    # now convert back to (l, d) for passing to function in PeelCoordinates module
    l_modified = dl_modified
    d_modified = d_modified + dl_modified

    return l_modified, d_modified


METHOD_NAME_TO_FUNCTION_POINT_CODE_TO_XYZ = {
    # "edge bisection": get_xyz_from_point_code_using_ancestry,
    # "uncorrected plane gridding": pe.get_xyz_from_point_code_using_peel_coordinates,
    "corrected plane gridding": get_xyz_from_point_code_using_corrected_plane_gridding,
    # "arc gridding": lambda pc: NotImplemented,
}


def get_stats_about_point_placements(pc_to_xyz):
    min_ds = []
    max_ds = []
    for pc, xyz in sorted(pc_to_xyz.items()):
        adj = get_adjacency_from_point_code(pc)
        ds = []
        for pc1 in adj.values():
            if pc1 is not None:
                xyz1 = pc_to_xyz.get(pc1)
                if xyz1 is not None:
                    d = mcm.xyz_distance(xyz, xyz1)
                    ds.append(d)
        min_ds.append(min(ds) if len(ds) > 0 else np.nan)
        max_ds.append(max(ds) if len(ds) > 0 else np.nan)
    plt.hist(min_ds, bins=50, label="mins", alpha=0.5)
    plt.hist(max_ds, bins=50, label="maxs", alpha=0.5)
    plt.legend()
    plt.show()



if __name__ == "__main__":
    pcs = ["A", "B"]
    for spc in ["C", "D", "E", "F"]:
        pcs += get_all_point_codes_from_ancestor_at_iteration(ancestor_pc=spc, iterations=6)

    # pcs = get_all_point_codes_at_iteration(iterations=5)

    for method_name, func_pc_to_xyz in METHOD_NAME_TO_FUNCTION_POINT_CODE_TO_XYZ.items():
        print(f"getting xyz coordinates using {method_name!r} method")
        xyzs = []
        pc_to_xyz = {}
        for pc in pcs:
            xyz = func_pc_to_xyz(pc)
            # print(f"{pc = }, {xyz = }")
            xyzs.append(xyz)
            pc_to_xyz[pc] = xyz

        # get_stats_about_point_placements(pc_to_xyz)

        # TODO sanity checks for a given method:
        # - the point's neighbors by adjacency must be the nearest neighbors
        # TODO check if any of the methods creates the same results

        ppl.plot_xyzs_on_sphere_3d(xyzs, labels=None)

        # plot a path
        pc_init, pc_final, direction = "C101010101", "C202020202", 3
        path_pcs = get_point_path(pc_init, pc_final, direction)
        ppl.plot_point_codes_on_half_peel_face_planes(path_pcs, face_name="CAKX", func_pc_to_xyz=func_pc_to_xyz, with_labels=False)
        ls, ds = pe.get_peel_coordinates_of_point_codes_on_face(path_pcs, face_name="CAKX", func_pc_to_xyz=func_pc_to_xyz)
        distances, angles = get_stepwise_path_distances_and_angles_2d(ls, ds)
        plot_distances_and_angles_2d(distances, angles)

        path_xyzs = [func_pc_to_xyz(pc) for pc in path_pcs]
        xs, ys, zs = zip(*path_xyzs)
        distances, angles_xy, angles_xz, angles_yz = get_stepwise_path_distances_and_angles_3d(xs, ys, zs)
        plot_distances_and_angles_3d(distances, angles_xy, angles_xz, angles_yz)

        print()
        input("press enter to continue")
        print()
    print("done")
