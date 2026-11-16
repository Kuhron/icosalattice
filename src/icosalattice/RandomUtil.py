import random
import numpy as np

from icosalattice.PointRepresentationAsFloat import point_float_to_code


def get_pseudorandom_point_codes(n_pcs, seed):
    random.seed(seed)
    ints = [random.randint(2, 11) for i in range(n_pcs)]
    floats = [random.random() for i in range(n_pcs)]
    xs = [a + b for a, b in zip(ints, floats)]
    pcs = [point_float_to_code(x, max_iterations=9, allow_clipping=True) for x in xs]
    return pcs

    # old stuff where I was trying to mess with the numbers myself

    # rs = [random.random() for i in range(n_pcs)]
    # xs = [random.randint(2, 120) for i in range(n_pcs)]
    # a = random.randint(13, 131)
    # b = random.randint(17, 41)

    # pcs = []
    # for i in range(n_pcs):
    #     j = a + b * i
    #     x = xs[i]
    #     r = rs[i]

    #     z = 1 + ((1+r) ** 234.919) % 77

    #     d = (x**(np.pi) % (1+x)**z + j**2 + np.pi) % 1
    #     k = int(x**3 / (j * 7 + 9.22023918) - (d*4)**2) % 10 + 2
    #     f = k + d
    #     pc = point_float_to_code(f, max_iterations=7, allow_clipping=True)
    #     pcs.append(pc)
    # return pcs


def get_test_floats_01(n):
    ns = np.pi * np.arange(1, n+1)
    xs = (ns**1.3 + ns**2.1) % 1
    return xs
