import pytest

import math
import numpy as np

from icosalattice.MathUtil import zigzag, zigzag_inverse


def test_zigzag_function():
    vs = [1/10, 1/3, 1/2, 1, 4/3, 3/2, 2, np.pi/4]
    xs = set()
    for a in vs:
        xs |= set(np.arange(0, 8.01, a))
    xs = sorted(xs)

    for a in vs:
        for x in xs:
            z = zigzag(x, a=a)
            n = math.floor(x/a)
            zi = zigzag_inverse(z, a=a, n=n)
            assert np.isclose(zi, x, rtol=1e-9), f"{a = :f}, {x = :f}, {z = :f}, {n = :f}, {zi = :f}"
