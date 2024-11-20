import icosalattice.Edges as ed
import icosalattice.StartingPoints as sp


def get_test_point_codes():
    edge_midpoints = ed.get_edge_midpoints()

    vertex_pcs = sp.STARTING_POINT_CODES
    edge_pcs = []
    for edge_name in edge_midpoints.keys():
        spc, direction = ed.get_ancestor_starting_point_and_direction_of_edge(edge_name)
        edge_pc = spc + direction
        edge_pcs.append(edge_pc)
    edge_pcs += ["C10101", "D3333", "E000303", "F22", "G111011", "H330333", "I202202", "J3003", "K1000001", "L22220202"]
    face_pcs = ["C10201", "D3320222", "E12", "F2023", "G11121", "H000102", "I01110113", "J3323", "K102", "L32121301"]

    pcs = vertex_pcs + edge_pcs + face_pcs
    return pcs


TEST_POINT_CODES = get_test_point_codes()
