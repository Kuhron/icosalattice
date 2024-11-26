import icosalattice.PeelCoordinates as pe
import icosalattice.FacePlaneDistortion as distort


def get_xyz_from_point_code_using_corrected_plane_gridding(pc, as_array=True):
    # modify the l and d coordinates to try to offset the distortion introduced when projecting from face plane back onto sphere surface
    spc, l_raw, d_raw = pe.get_peel_coordinates_from_point_code(pc)
    if l_raw > d_raw:
        # on upward-pointing face
        l_modified, d_modified = adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_upward_face(l_raw, d_raw)
    elif d_raw > l_raw:
        # on downward-pointing face
        l_modified, d_modified = adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_downward_face(l_raw, d_raw)
    else:
        # on direction-2 edge, check both faces
        l_modified_0, d_modified_0 = adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_upward_face(l_raw, d_raw)
        l_modified_1, d_modified_1 = adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_downward_face(l_raw, d_raw)
        assert abs(l_modified_0 - l_modified_1) < 1e-9, f"{l_modified_0} != {l_modified_1}"
        assert abs(d_modified_0 - d_modified_1) < 1e-9, f"{d_modified_0} != {d_modified_1}"
        l_modified = (l_modified_0 + l_modified_1) / 2
        d_modified = (d_modified_0 + d_modified_1) / 2

    return pe.get_xyz_from_peel_coordinates(spc, l_modified, d_modified, as_array=as_array)


def adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_upward_face(l_raw, d_raw):
    # D vector is not in the face plane! so need in terms of L and DL

    # convert (l, d) to (l, dl) for distorting
    new_dl_raw = d_raw
    new_l_raw = l_raw - d_raw
    dl_modified = distort.get_lp_proportion_from_theta_proportion(new_dl_raw)
    l_modified = distort.get_lp_proportion_from_theta_proportion(new_l_raw)

    # now convert back to (l, d) for passing to function in PeelCoordinates module
    d_modified = dl_modified
    l_modified = l_modified + dl_modified

    return l_modified, d_modified


def adjust_l_and_d_coordinates_for_corrected_plane_gridding_on_downward_face(l_raw, d_raw):
    # L vector is not in the face plane! so need in terms of DL and D
    
    # convert (l, d) to (dl, d) for distorting
    new_dl_raw = l_raw
    new_d_raw = d_raw - l_raw
    dl_modified = distort.get_lp_proportion_from_theta_proportion(new_dl_raw)
    d_modified = distort.get_lp_proportion_from_theta_proportion(new_d_raw)

    # now convert back to (l, d) for passing to function in PeelCoordinates module
    l_modified = dl_modified
    d_modified = d_modified + dl_modified

    return l_modified, d_modified
