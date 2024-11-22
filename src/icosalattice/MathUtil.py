# generic math functions that don't belong better somewhere else


import numpy as np


def round_off_unwanted_float_precision(x):
    # if x has a huge gap in order of magnitude between its second-least-significant digit and its least-significant
    # then the last one is likely float crap
    eps = 1e-9

    if x == 0:
        return x
    if x % eps == 0:
        return x
    
    a = 1
    while x < 1:
        a *= 2
        x *= 2
    
    rem = x % eps
    rem_over_eps = rem / eps
    # print(x, rem_over_eps)
    if rem_over_eps > 1 - 1e-4:
        # we need to ADD a small amount, round x up
        x2 = round(x, 9)
        assert x2 >= x
    elif rem_over_eps < 1e-4:
        x2 = round(x, 9)
        assert x2 <= x
    else:
        # keep the precision
        x2 = x
    return x2/a  # rescale since we multiplied it up


def get_vector_decomposition_coefficients(v, v1, v2):
    # potential optimizations, if needed:
    # - rewrite the matrix equation code as more direct expressions for a1/a2/a3 rather than inverting a matrix (just write the inverse yourself)
    # - - so that the math is done more directly with basic arithmetic operations in raw Python rather than NumPy

    # v = a1*v1 + a2*v2, solve for a1 and a2
    x1, y1, z1 = v1
    x2, y2, z2 = v2
    xp, yp, zp = v
    # solve equation for a1 and a2: [[xp] [yp]] = [[x1 x2] [y1 y2]] [[a1] [a2]]
    A = np.array([[x1, x2], [y1, y2]])
    Ainv = np.linalg.inv(A)
    a1, a2 = Ainv @ np.array([xp, yp])
    # verify solution works with z coordinates
    diff = zp - (a1*z1 + a2*z2)
    assert abs(diff) < 1e-9, "z coordinate verification of vector decomposition failed"

    # floats
    if -1e-9 < a1 < 0:
        a1 = 0
    if -1e-9 < a2 < 0:
        a2 = 0

    if not (0 <= a1 <= 1 and 0 <= a2 <= 1):
        raise InvalidVectorDecompositionException(f"vector decomposition should have coefficients between 0 and 1\ngot:\n  {v}\n= {a1} * {v1}\n+ {a2} * {v2}")
    a1 = float(a1)
    a2 = float(a2)
    return a1, a2


class InvalidVectorDecompositionException(Exception): pass
