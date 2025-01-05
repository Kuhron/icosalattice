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
from icosalattice.CoordinatesByPlaneGridding import get_xyz_from_point_code_using_corrected_plane_gridding
from icosalattice.CoordinatesByRThetaAdjustment import get_xyz_from_point_code_using_r_theta_adjustment
from icosalattice.GeneratePointCodes import get_all_point_codes_from_ancestor_at_iteration, get_all_point_codes_at_iteration
import icosalattice.PeelCoordinates as pe
import icosalattice.FacePlaneDistortion as distort
import icosalattice.PlotPointLocations as ppl
from icosalattice.Adjacency import get_adjacency_from_point_code
import icosalattice.MapCoordinateMath as mcm
from icosalattice.PlotPaths import plot_distances_and_angles_2d, plot_distances_and_angles_3d
from icosalattice.PointPaths import get_point_path, get_stepwise_path_distances_and_angles_2d, get_stepwise_path_distances_and_angles_3d



METHOD_NAME_TO_FUNCTION_POINT_CODE_TO_XYZ = {
    "ebs1": get_xyz_from_point_code_using_ancestry,  # "edge bisection"
    "upg1": pe.get_xyz_from_point_code_using_peel_coordinates,  # "uncorrected plane gridding"
    "cpg1": get_xyz_from_point_code_using_corrected_plane_gridding,  # "corrected plane gridding"
    "rta1": get_xyz_from_point_code_using_r_theta_adjustment,  # "r-theta adjustment"
}
CHOSEN_METHOD = "rta1"


def get_xyz_from_point_code(pc, as_array=True):
    f = METHOD_NAME_TO_FUNCTION_POINT_CODE_TO_XYZ[CHOSEN_METHOD]
    return f(pc, as_array=as_array)


def get_latlon_from_point_code(pc, as_array=True):
    xyz = get_xyz_from_point_code(pc)
    latlon = mcm.unit_vector_cartesian_to_latlon(*xyz, as_array=as_array)
    return latlon


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
    for spc in ["C", "E", "G", "I", "K"]:
        pcs += get_all_point_codes_from_ancestor_at_iteration(ancestor_pc=spc, iterations=6)

    # pcs = get_all_point_codes_at_iteration(iterations=5)

    method_name = "corrected plane gridding"
    func_pc_to_xyz = METHOD_NAME_TO_FUNCTION_POINT_CODE_TO_XYZ[method_name]
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
    # pc_init, pc_final, direction = "C101010101", "C202020202", 3
    pc_init, pc_final, direction = "C1110000", "C2220000", 3
    path_pcs = get_point_path(pc_init, pc_final, direction)
    ppl.plot_point_codes_on_half_peel_face_planes(path_pcs, face_name="CAKX", func_pc_to_xyz=func_pc_to_xyz, with_labels=False)
    ls, ds = pe.get_adjusted_peel_coordinates_of_point_codes_on_face(path_pcs, face_name="CAKX", func_pc_to_xyz=func_pc_to_xyz)
    distances, angles = get_stepwise_path_distances_and_angles_2d(ls, ds)
    plot_distances_and_angles_2d(distances, angles)

    path_xyzs = [func_pc_to_xyz(pc) for pc in path_pcs]
    xs, ys, zs = zip(*path_xyzs)
    distances, angles_xy, angles_xz, angles_yz = get_stepwise_path_distances_and_angles_3d(xs, ys, zs)
    plot_distances_and_angles_3d(distances, angles_xy, angles_xz, angles_yz)

