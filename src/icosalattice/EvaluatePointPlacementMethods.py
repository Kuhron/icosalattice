from functools import reduce
import numpy as np
import matplotlib.pyplot as plt

import icosalattice.GeneratePointCodes as gpc
import icosalattice.PlotPointLocations as ppl
from icosalattice.PointPaths import get_point_path, get_stepwise_path_distances_and_angles_2d
from icosalattice.PlotPaths import plot_distances_and_angles_2d
import icosalattice.PeelCoordinates as pe
from icosalattice.Adjacency import get_neighbors_of_point_code
from icosalattice.AnglesOnSphere import measure_angles_on_sphere
from icosalattice.DistancesOnSphere import distance_great_circle



def plot_point_distribution_on_one_face(func_pc_to_xyz, iterations=4):
    # observe how the points lie along lines or curves on the face
    # pcs = gpc.get_all_point_codes_from_ancestor_at_iteration(ancestor_pc="C", iterations=6)
    pcs = gpc.get_all_point_codes_on_face_at_iteration(face_name="CAKX", iterations=iterations, with_edges=True, with_trailing_zeros=True)
    ppl.plot_point_codes_on_half_peel_face_planes(pcs, face_name="CAKX", func_pc_to_xyz=func_pc_to_xyz, with_labels=False)
    # for future self searching for code that made plots:
    # CAKX_i4_{method}.png


def plot_point_distribution_on_northern_half_peels(func_pc_to_xyz, iterations=6):
    pcs = reduce((lambda x,y: x+y), [gpc.get_all_point_codes_from_ancestor_at_iteration(ancestor_pc=apc, iterations=iterations) for apc in "CEGIK"]) + ["A"]
    ppl.plot_point_codes_on_sphere_3d(pcs, func_pc_to_xyz=func_pc_to_xyz, with_labels=False)


def plot_trajectory_between_two_points_on_face_plane(pc_init, pc_final, direction, func_pc_to_xyz):
    pcs = get_point_path(pc_init, pc_final, direction)[:-1]
    ppl.plot_point_codes_on_half_peel_face_planes(pcs, face_name="CAKX", func_pc_to_xyz=func_pc_to_xyz, with_labels=False)
    ls, ds = pe.get_adjusted_peel_coordinates_of_point_codes_on_face(pcs, face_name="CAKX", func_pc_to_xyz=func_pc_to_xyz)
    distances, angles = get_stepwise_path_distances_and_angles_2d(ls, ds)
    plot_distances_and_angles_2d(distances, angles)
    # for future self searching for code that made plots:
    # C010100000_to_K010100000_steps_{method}.png
    # {pc_init}_to_{pc_final}_steps_{method}.png


def report_neighbor_angle_and_distance_statistics(func_pc_to_xyz, iterations=6):
    pcs = gpc.get_all_point_codes_at_iteration(iterations=iterations, with_trailing_zeros=True)
    pc_to_xyz = {pc: func_pc_to_xyz(pc) for pc in pcs}
    angles = []

    # debug: finding aberrantly large angles
    angles_larger_than_72 = {}

    distances = []
    for i, pc in enumerate(pcs):
        if i % 1000 == 0:
            print(f"getting distances and angles: {i = } / {len(pcs)}")
        xyz0 = pc_to_xyz[pc]
        neighbors = get_neighbors_of_point_code(pc)
        # print(pc, neighbors)
        these_angles_deg = []
        these_distances = []

        # there is double counting and redundant computation going on here
        # (since the distances are used in computing the angles, we compute the same angles more than once, etc.) 
        # but I don't care right now unless it takes very long to finish

        for i in range(len(neighbors)):
            n1 = neighbors[i]
            n2 = neighbors[(i+1) % len(neighbors)]
            xyz1 = pc_to_xyz[n1]
            xyz2 = pc_to_xyz[n2]
            angle1, angle0, angle2 = measure_angles_on_sphere(xyz1, xyz0, xyz2)
            angle0_deg = angle0 * 180/np.pi
            
            # debug: finding aberrantly large angles
            if angle0_deg > 72 + 1e-9:
                a = round(angle0_deg, 6)
                if a not in angles_larger_than_72:
                    angles_larger_than_72[a] = []
                angles_larger_than_72[a].append((n1, pc, n2))

            distance = distance_great_circle(xyz0, xyz1)
            these_angles_deg.append(angle0_deg)
            these_distances.append(distance)
        angles += these_angles_deg
        distances += these_distances
    
    print("\nlarge angles:\n" + ("\n".join(f"  {a:f} deg: {tup}" for a, tups in sorted(angles_larger_than_72.items(), reverse=True) for tup in tups) if len(angles_larger_than_72) > 0 else "(none)") + "\n")
    
    print(f"{min(angles) = :f}")
    print(f"mean(angles) = {np.mean(angles):f}")
    print(f"median(angles) = {np.median(angles):f}")
    print(f"{max(angles) = :f}\n")

    print(f"{min(distances) = :f}")
    print(f"mean(distances) = {np.mean(distances):f}")
    print(f"median(distances) = {np.median(distances):f}")
    print(f"{max(distances) = :f}\n")

    plt.hist(angles, bins=100)
    plt.title("angles")
    plt.show()
    # for future self searching for code that made plots:
    # neighbor_angle_distribution_{method}.png

    plt.hist(distances, bins=100)
    plt.title("distances")
    plt.show()
    # for future self searching for code that made plots:
    # neighbor_distance_distribution_{method}.png
