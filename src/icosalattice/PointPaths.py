# functions about paths between two points
# or series of points forming a path


import icosalattice.PointCodeArithmetic as pca


def get_point_path(pc_init, pc_final, direction):
    pcs = [pc_init]
    last_pc_added = pc_init
    pcs_seen = {pc_init}
    while last_pc_added != pc_final:
        next_pc = pca.add_direction_to_point_code(last_pc_added, direction)
        if next_pc in pcs_seen:
            raise RuntimeError(f"loop detected: point code {next_pc} is already in the path in the {direction} direction from {pc_init}, so the destination {pc_final} will never be reached")
        pcs.append(next_pc)
        pcs_seen.add(next_pc)
    return pcs

