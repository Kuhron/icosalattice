import math

import icosalattice.MapCoordinateMath as mcm
from icosalattice.UnitSpherePoint import UnitSpherePoint
from icosalattice.ConstantMakerDecorator import constant_maker


# constants that can be defined before functions
MID_LAT_DEG = math.atan(1/2) * 180/math.pi
STARTING_POINT_CODES = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"]


def get_starting_points_latlon_named():
    icosahedron_original_points_latlon = {
        # north pole
        "NP": (90, 0),
        # north ring of five points, star with a point at lon 0
        "NR0": (MID_LAT_DEG, 0), "NRp72": (MID_LAT_DEG, 72), "NRp144": (MID_LAT_DEG, 144), "NRm72": (MID_LAT_DEG, -72), "NRm144": (MID_LAT_DEG, -144),
        # south ring of five points, star with a point at lon 180
        "SR180": (-MID_LAT_DEG, 180), "SRp108": (-MID_LAT_DEG, 108), "SRp36": (-MID_LAT_DEG, 36), "SRm108": (-MID_LAT_DEG, -108), "SRm36": (-MID_LAT_DEG, -36),
        # south pole
        "SP": (-90, 0),
    }
    return icosahedron_original_points_latlon


def get_starting_points_adjacency_named():
    original_points_neighbors_by_name = {
        # start with "north" neighbor, i.e., directly left on peel-rectangle representation (for poles the order should still obey counterclockwise, but starting point doesn't matter)
        # poles
        "NP": ["NR0", "NRp72", "NRp144", "NRm144", "NRm72"],  # going eastward (counterclockwise)
        "SP": ["SR180", "SRp108", "SRp36", "SRm36", "SRm108"],  # going westward (counterclockwise)

        # peel 0
        "NR0": ["NP", "NRm72", "SRm36", "SRp36", "NRp72"],
        "SRp36": ["NR0", "SRm36", "SP", "SRp108", "NRp72"],

        # peel 1
        "NRp72": ["NP", "NR0", "SRp36", "SRp108", "NRp144"],
        "SRp108": ["NRp72", "SRp36", "SP", "SR180", "NRp144"],

        # peel 2
        "NRp144": ["NP", "NRp72", "SRp108", "SR180", "NRm144"],
        "SR180": ["NRp144", "SRp108", "SP", "SRm108", "NRm144"],

        # peel 3
        "NRm144": ["NP", "NRp144", "SR180", "SRm108", "NRm72"],
        "SRm108": ["NRm144", "SR180", "SP", "SRm36", "NRm72"],

        # peel 4
        "NRm72": ["NP", "NRm144", "SRm108", "SRm36", "NR0"],
        "SRm36": ["NRm72", "SRm108", "SP", "SRp36", "NR0"],
    }
    assert len(original_points_neighbors_by_name) == 12 and all(len(vals) == 5 for vals in original_points_neighbors_by_name.values())
    # check transitivity of neighborliness, since I input the lists manually
    for point_name, neighbors in original_points_neighbors_by_name.items():
        for neigh in neighbors:
            assert point_name in original_points_neighbors_by_name[neigh], "intransitive adjacency with {} and {}".format(point_name, neigh)

    return original_points_neighbors_by_name


def get_original_points_order_by_name():
    original_points_order_by_name = [
        "NP", "SP",  # poles
        "NR0", "SRp36",  # peel 0
        "NRp72", "SRp108",  # peel 1
        "NRp144", "SR180",  # peel 2
        "NRm144", "SRm108",  # peel 3
        "NRm72", "SRm36",  # peel 4
    ]
    return original_points_order_by_name


def get_starting_points():
    # print("getting starting icosa points")
    icosahedron_original_points_latlon = get_starting_points_latlon_named()
    original_points_neighbors_by_name = get_starting_points_adjacency_named()

    # put them in this ordering convention:
    # north pole first, south pole second, omit these from all expansion operations, by only operating on points[2:] (non-pole points)
    # new points are appended to the point list in the order they are created
    # order the neighbors in the following order, and only bisect the edges from each point to the first three:
    # - [north, west, southwest, (others)]. order of others doesn't matter that much, can just keep going counterclockwise
    # ordering of neighbors for poles thus doesn't matter as that list will never be used for expansion

    original_points_order_by_name = get_original_points_order_by_name()

    # keep the point objects in a single array that can be indexed by point index
    # the rest of the data, i.e., the adjacencies dictionary, should be all in terms of integer indices that refer to the points array

    ordered_points = []
    adjacencies_by_point_index = [None for i in range(12)]

    # place original points in the list
    for point_number, p_name in enumerate(original_points_order_by_name):
        point_index = len(ordered_points)
        p_latlon = icosahedron_original_points_latlon[p_name]
        p_xyz = mcm.unit_vector_latlon_to_cartesian(*p_latlon)
        coords_dict = {"xyz": p_xyz, "latlondeg": p_latlon}
        usp = UnitSpherePoint(coords_dict, point_number)
        ordered_points.append(usp)
    assert len(ordered_points) == 12, "initial icosa needs 12 vertices"

    # add their neighbors by index
    for point_index in range(len(ordered_points)):
        point_name = original_points_order_by_name[point_index]
        neighbor_names = original_points_neighbors_by_name[point_name]
        neighbor_indices = [original_points_order_by_name.index(name) for name in neighbor_names]
        adjacencies_by_point_index[point_index] = neighbor_indices
        # print("adjacencies now:\n{}\n".format(adjacencies_by_point_index))

    # print("-- done getting initial icosa points")
    return ordered_points, adjacencies_by_point_index


@constant_maker("STARTING_POINTS")
def get_starting_points_immutable():
    ordered_points, adj = get_starting_points()
    assert type(ordered_points) is list
    assert all(type(x) is UnitSpherePoint for x in ordered_points)
    ordered_points = tuple(ordered_points)
    assert type(adj) is list
    new_adj = ()
    for x in adj:
        assert type(x) is list
        assert all(type(y) is int for y in x)
        x_tup = tuple(x)
        extend_tup = (x_tup,)
        new_adj = new_adj + extend_tup
    return (ordered_points, new_adj)


def get_starting_point_neighbor_identity(point_number):
    # for 0 and 1 (the poles) this is still weird, it's not clear what the directions (L,DL,D,R,UR,U) would mean for them, ill-defined like east of the south pole
    # but for the other 10 starting points, there are five neighbors but one of them acts like two directions
    # e.g. on the northern ring, from the perspective of the peel below (west of) the point, the L neighbor is the north pole
    # but from the perspective of the peel above (east of) the point, the U neighbor is the north pole
    d = {}
    assert type(point_number) is int, point_number
    assert 2 <= point_number < 12, "invalid point for neighbor identity: {}".format(point_number)
    ring = get_starting_point_ring(point_number)
    if ring == "northern_ring":
        return ("L", "U")
    elif ring == "southern_ring":
        return ("D", "R")
    else:
        raise ValueError("invalid ring {}".format(ring))


def get_starting_point_ring(starting_point):
    original_points_order_by_name = get_original_points_order_by_name()
    ring_code = original_points_order_by_name[starting_point][:2]
    if ring_code == "NP":
        return "north_pole"
    elif ring_code == "SP":
        return "south_pole"
    elif ring_code == "NR":
        return "northern_ring"
    elif ring_code == "SR":
        return "southern_ring"
    else:
        raise ValueError("invalid ring code {}".format(ring_code))


@constant_maker("STARTING_DIRECTIONAL_DICT")
def get_starting_point_code_directional_dict():
    # this will help with finding the directional parent for a given point code
    # A and B don't have directions (L,DL,D = 1,2,3) coming from them
    return {
        "C": {"1": "A", "2": "K", "3": "L"},
        "D": {"1": "C", "2": "L", "3": "B"},
        "E": {"1": "A", "2": "C", "3": "D"},
        "F": {"1": "E", "2": "D", "3": "B"},
        "G": {"1": "A", "2": "E", "3": "F"},
        "H": {"1": "G", "2": "F", "3": "B"},
        "I": {"1": "A", "2": "G", "3": "H"},
        "J": {"1": "I", "2": "H", "3": "B"},
        "K": {"1": "A", "2": "I", "3": "J"},
        "L": {"1": "K", "2": "J", "3": "B"},
    }


# constants that need to be defined after the function that creates them
STARTING_POINTS_AND_ADJACENCY = get_starting_points_immutable(calling_to_create_constant=True)  # since this is called way too many times otherwise, just initialize it as a global constant that can be accessed by further functions, e.g. base case for recursive adjacency algorithm
STARTING_POINTS, STARTING_ADJACENCY = STARTING_POINTS_AND_ADJACENCY
STARTING_DIRECTIONAL_DICT = get_starting_point_code_directional_dict(calling_to_create_constant=True)
