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
        "AC": (high_mid, 0),
        "AE": (high_mid, 72),
        "AG": (high_mid, 144),
        "AI": (high_mid, -144),
        "AK": (high_mid, -72),
        "BD": (low_mid, 36),
        "BF": (low_mid, 108),
        "BH": (low_mid, 180),
        "BJ": (low_mid, -108),
        "BL": (low_mid, -36),
        "CE": (high_mid_arc, 36),
        "EG": (high_mid_arc, 108),
        "GI": (high_mid_arc, 180),
        "IK": (high_mid_arc, -108),
        "KC": (high_mid_arc, -36),
        "DF": (low_mid_arc, 72),
        "FH": (low_mid_arc, 144),
        "HJ": (low_mid_arc, -144),
        "JL": (low_mid_arc, -72),
        "LD": (low_mid_arc, 0),
        "CD": (0, 18),
        "DE": (0, 54),
        "EF": (0, 90),
        "FG": (0, 126),
        "GH": (0, 162),
        "HI": (0, -162),
        "IJ": (0, -126),
        "JK": (0, -90),
        "KL": (0, -54),
        "LC": (0, -18),
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
