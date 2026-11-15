import pytest
import numpy as np

import icosalattice.PeelCoordinates as pe
from icosalattice.TriangularPeelCoordinates import get_ld_from_ack, get_ack_from_ld

from TestUtil import TEST_POINT_CODES



def test_ld_ack_conversion():
    for pc in TEST_POINT_CODES:
        sld_raw = pe.get_raw_peel_coordinates_from_point_code(pc)
        spc, l, d = sld_raw
        a, c, k = get_ack_from_ld(l, d)
        l2, d2 = get_ld_from_ack(a, c, k)
        a2, c2, k2 = get_ack_from_ld(l2, d2)
        assert np.isclose(l, l2, atol=1e-9)
        assert np.isclose(d, d2, atol=1e-9)
        assert np.isclose(a, a2, atol=1e-9)
        assert np.isclose(c, c2, atol=1e-9)
        assert np.isclose(k, k2, atol=1e-9)
