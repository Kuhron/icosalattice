import pytest

import icosalattice.StartingPoints as sp
import icosalattice.MapCoordinateMath as mcm


def test_starting_points_equidistant():
    points, adjacencies = sp.STARTING_POINTS_AND_ADJACENCY
    distances = set()
    for i, p in enumerate(points):
        neighbors = [points[j] for j in adjacencies[i]]
        for p2 in neighbors:
            d = mcm.xyz_distance(p.xyz(), p2.xyz())
            distances.add(d)
    assert max(distances) - min(distances) < 1e-9, f"starting points not equidistant, got min {min(distances)} and max {max(distances)}"
