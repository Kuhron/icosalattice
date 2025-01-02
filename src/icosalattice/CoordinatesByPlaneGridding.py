import icosalattice.PeelCoordinates as pe
import icosalattice.FacePlaneDistortion as distort
from icosalattice.TriangularPeelCoordinates import adjust_ld_using_lp_transformation_in_triangle_coordinates


def get_xyz_from_point_code_using_corrected_plane_gridding(pc, as_array=True):
    # modify the l and d coordinates to try to offset the distortion introduced when projecting from face plane back onto sphere surface
    spc, l_raw, d_raw = pe.get_raw_peel_coordinates_from_point_code(pc)
    l_modified, d_modified = adjust_ld_using_lp_transformation_in_triangle_coordinates(l_raw, d_raw)
    sld = spc, l_modified, d_modified
    return pe.get_xyz_from_adjusted_peel_coordinates(sld, as_array=as_array)
