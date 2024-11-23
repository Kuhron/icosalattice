# different methods for converting a point code into coordinates
# Warning: these different methods may generate different lattices!

# name of original method for locating child points: Edge Bisection

# TODO one idea for redoing where the child points are:
# - method name: Corrected Plane Gridding
# - use peel coordinates as primary
# - but dilate where the l and d proportions are by using the trig formula for lp (so the same delta l moves less near the center and more near the corners)
# - draw triangular grid on the face plane based solely on these lp points on the three edges, so the grid lines are straight on the face plane
# - project the resulting points to the sphere
# - see what they look like and how far apart they are, hopefully the lp adjustment makes them relatively uniformly spaced at least within a face
# - and hopefully the lines put the points in neat rows on the sphere surface, and most or all refraction occurs at actual face boundaries
# - (want no Sierpinski artifacts in the distribution of points within a face)

# TODO another idea for child point locations:
# - method name: Arc Gridding
# - divide the sphere edges into n equal segments based on which iteration you are at (e.g. 8 segments)
# - then connect corresponding points along these edges using great circle paths
# - hopefully the gc paths meet at sixfold vertices
# - and I suspect this method is equivalent to Grid Distortion

# the method of doing the gridding on the face plane uniformly and then projecting that out onto the sphere, introducing distortion,
# - could be called Uncorrected Plane Gridding


from icosalattice.CoordinatesByAncestry import get_xyz_from_point_code_using_ancestry
from icosalattice.GeneratePointCodes import get_all_point_codes_from_ancestor_at_iteration


METHOD_NAME_TO_FUNCTION_POINT_CODE_TO_XYZ = {
    "bisection": get_xyz_from_point_code_using_ancestry,
    # "uncorrected plane gridding": lambda pc: NotImplemented,
    # "corrected plane gridding": lambda pc: NotImplemented,
    # "arc gridding": lambda pc: NotImplemented,
}


if __name__ == "__main__":
    pcs = get_all_point_codes_from_ancestor_at_iteration(ancestor_pc="C01", iterations=6)
    print(pcs)

    for method_name, func in METHOD_NAME_TO_FUNCTION_POINT_CODE_TO_XYZ.items():
        print(f"getting xyz coordinates using {method_name!r} method")
        xyzs = []
        for pc in pcs:
            xyz = func(pc)
            print(f"{pc = }, {xyz = }")
            xyzs.append(xyz)
        # TODO run some stats on the points, like histogram of distances to neighbors, and plot them
        print()
        input("press enter to continue")
        print()
    print("done")
