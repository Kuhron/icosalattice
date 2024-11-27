# measure how much the lattice points are over/undersampled across the globe


import numpy as np
import matplotlib.pyplot as plt
# from scipy.spatial import distance_matrix

import icosalattice.GeneratePointCodes as gpc
from icosalattice.CoordinatesByPlaneGridding import get_xyz_from_point_code_using_corrected_plane_gridding
import icosalattice.MapCoordinateMath as mcm
from icosalattice.PlotDataOnMap import plot_variable_interpolated_from_dict
from icosalattice.Adjacency import get_adjacency_from_point_code


pcs = gpc.get_all_point_codes_at_iteration(iterations=6)
pc_to_xyz = {pc: get_xyz_from_point_code_using_corrected_plane_gridding(pc) for pc in pcs}
n = len(pcs)

average_distance_to_neighbors = {}
for pc in pcs:
    xyz = pc_to_xyz[pc]
    adj = get_adjacency_from_point_code(pc)
    neighbor_xyzs = [pc_to_xyz[pc1] for pc1 in adj.values() if pc1 is not None]  # can't get any data at the poles this way since they have all-None adjacency but whatever
    neighbor_distances = [mcm.xyz_distance(xyz, xyz1) for xyz1 in neighbor_xyzs]
    if len(neighbor_distances) > 0:
        average_distance_to_neighbors[pc] = np.mean(neighbor_distances)

plot_variable_interpolated_from_dict(average_distance_to_neighbors, dots_per_degree=1)
