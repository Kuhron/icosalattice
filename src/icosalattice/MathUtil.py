# generic math functions that don't belong better somewhere else

import math
import numpy as np



def round_off_unwanted_float_precision(x, epsilon=1e-9, min_abs=1e-12):
    # if x has a huge gap in order of magnitude between its second-least-significant digit and its least-significant
    # then the last one is likely float crap
    
    # x is zero
    if x == 0:
        return 0.0
    if abs(x) < min_abs:
        return 0.0
    
    # x has no float crap beyond the precision we want
    if x % epsilon == 0:
        return x
    
    neg = x < 0
    if neg:
        x = -x
    
    a = 1
    while x < 1:
        a *= 2
        x *= 2
    
    rem = x % epsilon
    rem_over_eps = rem / epsilon
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
    
    res = x2/a  # rescale since we multiplied it up
    if neg:
        return -res
    else:
        return res


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


def mod(x, y, map_zero_up=False):
    m = x % y
    if map_zero_up and m == 0:
        return y
    else:
        return m


def zigzag(x, a):
    # write logic only for a=1 case for simplicity and comprehensibility
    x /= a

    n = math.floor(x)
    n_odd = n % 2 == 1
    res = mod(x * (-1)**n, 1, map_zero_up=n_odd)

    res *= a
    return res


def zigzag_inverse(x, a, n):
    # write logic only for a=1 case for simplicity and comprehensibility
    x /= a

    res = x * (-1)**n + n + mod(n, 2)

    res *= a
    return res


class InvalidVectorDecompositionException(Exception): pass
