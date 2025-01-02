import icosalattice.Edges as ed
import icosalattice.StartingPoints as sp
from icosalattice.PointRepresentationAsFloat import point_float_to_code


def get_test_point_codes():
    edge_midpoints = ed.get_edge_midpoints()

    vertex_pcs = sp.STARTING_POINT_CODES
    edge_pcs = []
    for edge_name in edge_midpoints.keys():
        spc, direction = ed.get_ancestor_starting_point_and_direction_of_edge(edge_name)
        edge_pc = spc + direction
        edge_pcs.append(edge_pc)
    edge_pcs += ["C10101", "D3333", "E000303", "F22", "G111011", "H330333", "I202202", "J3003", "K1000001", "L2222022"]
    face_pcs = ["C10201", "D3320222", "E12", "F2023", "G11121", "H000102", "I0111013", "J3323", "K102", "L3212101"]

    pseudorandom_pcs = []
    for x in range(13, 13+29*18+1, 29):
        r = x**3.14159 % 1
        i = int(x**3 / 7.22023918 - (r*4)**2) % 10 + 2
        f = i + r
        pc = point_float_to_code(f, max_iterations=7, allow_clipping=True)
        pseudorandom_pcs.append(pc)
    
    various_pcs = ["H2223213", "K0010231", "J1203121", "D22023", "F1121", "G101021", "I3333202"]

    no_trailing_zero_pcs = vertex_pcs + edge_pcs + face_pcs + pseudorandom_pcs + various_pcs
    trailing_zero_pcs = [pc + "0" for pc in no_trailing_zero_pcs]

    pcs = no_trailing_zero_pcs + trailing_zero_pcs
    return pcs


TEST_POINT_CODES = get_test_point_codes()
