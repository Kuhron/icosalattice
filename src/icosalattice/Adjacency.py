import icosalattice.PointCodeArithmetic as pca
import icosalattice.StartingPoints as sp


def get_adjacency_from_point_code(pc, fix_edge_polarity=True):
    neighL = pca.add_direction_to_point_code(pc, 1, fix_edge_polarity=fix_edge_polarity)
    neighDL = pca.add_direction_to_point_code(pc, 2, fix_edge_polarity=fix_edge_polarity)
    neighD = pca.add_direction_to_point_code(pc, 3, fix_edge_polarity=fix_edge_polarity)
    neighR = pca.add_direction_to_point_code(pc, -1, fix_edge_polarity=fix_edge_polarity)
    neighUR = pca.add_direction_to_point_code(pc, -2, fix_edge_polarity=fix_edge_polarity)
    neighU = pca.add_direction_to_point_code(pc, -3, fix_edge_polarity=fix_edge_polarity)
    adj = {1: neighL, 2: neighDL, 3: neighD, -1: neighR, -2: neighUR, -3: neighU}
    # print(f"{pc=} has {adj=}")
    return adj


def get_neighbors_of_point_code(pc):
    # go in counterclockwise direction, which is how the [1, 2, 3, -1, -2, -3] vectors go for a point with 6 neighbors
    # around the north pole, this is the order CEGIK
    # around the south pole, this is the order LJHFD
    if pc[0] in sp.POLES:
        pca.validate_point_code(pc)
        iterations = len(pc) - 1
        if pc[0] == sp.NORTH_POLE:
            return [x + "1"*iterations for x in sp.NORTHERN_RING]
        elif pc[0] == sp.SOUTH_POLE:
            return [x + "3"*iterations for x in sp.SOUTHERN_RING[::-1]]
        else:
            raise Exception(f"invalid pole; shouldn't happen")
    else:
        adj = get_adjacency_from_point_code(pc)
        res = []
        for direction in [1, 2, 3, -1, -2, -3]:
            val = adj[direction]
            if val is not None:
                res.append(val)
        if all(x == "0" for x in pc[1:]):
            n_neighbors_expected = 5
        else:
            n_neighbors_expected = 6
        assert len(res) == n_neighbors_expected, f"{pc = } has wrong number of neighbors: {adj}"
        return res
