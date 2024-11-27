import icosalattice.PointCodeArithmetic as pca


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

