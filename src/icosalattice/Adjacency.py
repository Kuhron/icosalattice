import icosalattice.PointCodeArithmetic as pca


def get_adjacency_from_point_code(pc):
    neighL = pca.add_direction_to_point_code(pc, 1)
    neighDL = pca.add_direction_to_point_code(pc, 2)
    neighD = pca.add_direction_to_point_code(pc, 3)
    neighR = pca.add_direction_to_point_code(pc, -1)
    neighUR = pca.add_direction_to_point_code(pc, -2)
    neighU = pca.add_direction_to_point_code(pc, -3)
    adj = {1: neighL, 2: neighDL, 3: neighD, -1: neighR, -2: neighUR, -3: neighU}
    # print(f"{pc=} has {adj=}")
    return adj

