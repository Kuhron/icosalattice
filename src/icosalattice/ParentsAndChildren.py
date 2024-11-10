from icosalattice.Adjacency import get_adjacency_from_point_code
import icosalattice.BoxCornerMapping as bc
import icosalattice.IcosahedronMath as icm
import icosalattice.PointCodeArithmetic as arith


def get_parent_from_point_code(pc):
    if len(pc) == 1:
        return None
    if pc[-1] == "0":
        return pc[:-1]  # treat it as its own parent when it stays in the same place
    return arith.strip_trailing_zeros(pc[:-1])


def get_directional_parent_from_point_code(pc):
    if pc[-1] == "0":
        return pc[:-1]  # treat it as its own dpar when it stays in the same place
    if len(pc) == 1:
        return None
    
    reference_peel = bc.get_peel_containing_point_code(pc)

    if bc.point_code_is_in_reversed_polarity_encoding(pc):
        raise ValueError(f"got reverse-encoded {pc=}, please correct its encoding somewhere that we know which peel's perspective we are in")
    
    # faster than using adjacency and getting dpar at same index as this point's child_index
    # convert it to the C-D peel and then convert the answer back to the correct peel
    pc, peel_offset = arith.normalize_peel(pc)
    assert pc[0] in ["C", "D"], "peel normalization failed"

    # cd_dpar = bc.get_directional_parent_from_point_code_using_box_corner_mapping(pc)
    birth_direction = get_birth_direction_from_point_code(pc)
    par = get_parent_from_point_code(pc)
    cd_dpar = arith.add_direction_to_point_code(par, birth_direction)
    
    # keep trailing zeros during the recursive calls to box corner mapping
    # but no longer need them here now that we have our final answer
    cd_dpar = arith.strip_trailing_zeros(cd_dpar)

    # similarly undo reversed-polarity encoding if we got one of those
    if bc.point_code_is_in_reversed_polarity_encoding(cd_dpar):
        rev_cd_dpar = bc.correct_reversed_edge_polarity(cd_dpar, reference_peel)
        # print(f"reversed {cd_dpar} to {rev_cd_dpar}")
        cd_dpar = rev_cd_dpar

    dpar = arith.apply_peel_offset(cd_dpar, peel_offset)
    # print(f"offset {cd_dpar} to {dpar}")
    return dpar


def get_child_from_point_code(pc, child_index, iteration):
    if pc in ["A", "B"]:
        raise ValueError(f"point {pc} cannot have children")
    verify_can_have_children_from_point_code(pc, iteration)
    # pad with zeros if needed to get child at this iteration
    pc0 = pc.ljust(iteration, "0")
    child_index_str = str(child_index + 1)  # "0" is reserved for not moving and then making later child
    assert child_index_str in ["1", "2", "3"], child_index_str
    return pc0 + child_index_str


def get_children_from_point_code(pc):
    # including the null child where we keep the point in the same place
    if pc[0] in ["A", "B"]:
        # pole can't have children
        if all(y == "0" for y in pc[1:]):
            # it's actually the pole
            # but we need to return the pole itself at this iteration (pseudo-child)
            return [pc + "0"]
        else:
            raise ValueError(f"got reverse-encoded or invalid point code {pc}")
    else:
        return [pc + x for x in "0123"]


def get_birth_direction_from_point_code(pc):
    return int(pc[-1])


def verify_can_have_children_from_point_code(pc, iteration):
    # iteration must be greater than point's born iteration
    if pc in ["A", "B"]:
        raise ValueError(f"point {pc} cannot have children")
    return iteration > icm.get_iteration_born_from_point_code(pc)

