import numpy as np

import icosalattice.StartingPoints as sp
from icosalattice.UnitSpherePoint import UnitSpherePoint


MID_ARC_LAT = np.arctan((1+5**0.5)/2 -1) * 180/np.pi  # where the edge CE peaks in latitude


def get_edge_midpoints():
    ring_lat = sp.MID_LAT_DEG
    high_mid_arc = MID_ARC_LAT
    high_mid = 1/2 * (90 + ring_lat)
    low_mid = -high_mid
    low_mid_arc = -MID_ARC_LAT

    midpoints_latlon = {
        "CA": (high_mid, 0),
        "EA": (high_mid, 72),
        "GA": (high_mid, 144),
        "IA": (high_mid, -144),
        "KA": (high_mid, -72),
        "DB": (low_mid, 36),
        "FB": (low_mid, 108),
        "HB": (low_mid, 180),
        "JB": (low_mid, -108),
        "LB": (low_mid, -36),
        "EC": (high_mid_arc, 36),
        "GE": (high_mid_arc, 108),
        "IG": (high_mid_arc, 180),
        "KI": (high_mid_arc, -108),
        "CK": (high_mid_arc, -36),
        "FD": (low_mid_arc, 72),
        "HF": (low_mid_arc, 144),
        "JH": (low_mid_arc, -144),
        "LJ": (low_mid_arc, -72),
        "DL": (low_mid_arc, 0),
        "DC": (0, 18),
        "ED": (0, 54),
        "FE": (0, 90),
        "GF": (0, 126),
        "HG": (0, 162),
        "IH": (0, -162),
        "JI": (0, -126),
        "KJ": (0, -90),
        "LK": (0, -54),
        "CL": (0, -18),
    }
    return {k: UnitSpherePoint.from_latlondeg(*v) for k,v in midpoints_latlon.items()}


def get_ancestor_starting_point_and_direction_of_edge(endpoint_codes):
    d = sp.STARTING_DIRECTIONAL_DICT
    a,b = endpoint_codes
    ab = None
    ba = None
    if a in d:
        for direction, pc in d[a].items():
            if pc == b:
                ab = direction
    if b in d:
        for direction, pc in d[b].items():
            if pc == a:
                ba = direction

    assert (ab is not None) + (ba is not None) == 1, f"should have exactly one direction of edge {endpoint_codes} in starting dict"
    if ab is not None:
        spc = a
        direction = ab
    elif ba is not None:
        spc = b
        direction = ba
    else:
        raise Exception("impossible")
    return spc, direction
