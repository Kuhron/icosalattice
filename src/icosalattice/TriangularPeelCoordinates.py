import numpy as np
import random

import icosalattice.FacePlaneDistortion as distort


# conversion between seeing a half-peel as two equilateral triangles (triangle representation) vs two halves of a square (square representation)
# e.g. for faces CAKX and CXKL, origin is at C, l_hat vector is C to A, d_hat vector is C to L
# in each representation, length of l_hat and d_hat is 1
# in triangle representation, C = (0, 0), A = (-1/2, +R3/2), K = (-1, 0), L = (-1/2, -R3/2)
# in square representation, C = (0, 0), A = (-1, 0), K = (-1, -1), L = (0, -1)
# these matrices transform (x,y) coordinates (NOT (l,d)) from one of these ways of viewing the half-peel into the other
R3 = 3**0.5
MATRIX_TRIANGLE_TO_SQUARE = np.array([[1.0, -R3/3], [1.0, +R3/3]])
MATRIX_SQUARE_TO_TRIANGLE = np.array([[1/2, 1/2], [-R3/2, +R3/2]])

# I refer to MATRIX_TRIANGLE_TO_SQUARE also as A_t>q, and MATRIX_SQUARE_TO_TRIANGLE also as A_q>t

# on an upward-pointing triangular face plane, the distances a, c, k are defined as
# # the length along the edges from A, C, K (respectively) to the line that runs through the point (p)
# # and is parallel to the edge opposite A, C, K (respectively)

# diagram of triangle with these lines shows that, because the face plane is divided into equilateral triangles and rhombuses,
# # letting the segments on the edge from A to K be labeled x1, x2, x3,
# # then the segments from K to C are x2, x3, x1, and the segments from C to A are x3, x1, x2
# # so we have x1 + x2 + x3 = 1
# # and a + c + k = 2*(x1 + x2 + x3) = 2

# in the square representation:
# # l_q (the "left" coordinate of p) is the vector from C to the point on the edge CA which is directly above p
# # d_q (the "down" coordinate of p) is the vector from C to the point on the edge CL which is directly to the right of p

# in the triangle representation:
# # c = norm(A_q>t @ l_q)
# # k = 1 - norm(A_q>t @ d_q)
# # a = 2 - c - k

# on a downward-pointing triangular face plane, define a, c, k as negative
# # in parallel direction to a, c, k on the corresponding upward-pointing face plane
# # so here a is measured from L into the face, c is measured from K into the face, and k is measured from C into the face
# make sure the a,c,k are floats so we can have positive and negative zero which mean different things
POS_ZERO = +0.0
NEG_ZERO = -0.0

# in the square representation, l_q and d_q are defined the same way as for the upward-pointing face
# in the triangle representation:
# # |k| = norm(A_q>t @ d_q)
# # |c| = 1 - norm(A_q>t @ l_q)
# # |a| = 2 - |c| - |k|

# TODO try adjusting all of a, c, and k using the l_p(theta) transformation, and see if it
# - (1) works at all (maybe the lines don't all three meet at the same point anymore)
# - - 2024-11-26: the transformed a,c,k no longer obey a+c+k = 2 (it's always slightly less than 2, and too different to be due to float rounding)
# - - so right now it looks like this DOESN'T work
# - - maybe I can just scale them so that they add to 2 again? hack
# - (2) makes the points nicely spaced on the sphere (which is the point of doing all this peel coordinate distortion in the first place)


def convert_xy_triangle_to_square(xy):
    return MATRIX_TRIANGLE_TO_SQUARE @ xy


def convert_xy_square_to_triangle(xy):
    return MATRIX_SQUARE_TO_TRIANGLE @ xy


def get_ack_from_ld(l_coord, d_coord):
    if l_coord > d_coord:
        # upward face
        return get_ack_from_ld_on_upward_face(l_coord, d_coord)
    elif l_coord < d_coord:
        # downward face
        return get_ack_from_ld_on_downward_face(l_coord, d_coord)
    else:
        # edge between them
        a_up, c_up, k_up = get_ack_from_ld_on_upward_face(l_coord, d_coord)
        a_dn, c_dn, k_dn = get_ack_from_ld_on_downward_face(l_coord, d_coord)
        assert np.isclose(a_up, 1.0, atol=1e-9), a_up
        assert np.isclose(a_dn, -1.0, atol=1e-9), a_dn
        assert np.isclose(c_dn, -k_up, atol=1e-9), (c_dn, k_up)
        assert np.isclose(k_dn, -c_up, atol=1e-9), (k_dn, c_up)
        # print(f"on upward face  : a={a_up:.4f}, c={c_up:.4f}, k={k_up:.4f}")
        # print(f"on downward face: a={a_dn:.4f}, c={c_dn:.4f}, k={k_dn:.4f}")
        
        # just use the upward ones
        return a_up, c_up, k_up


def get_ack_from_ld_on_upward_face(l_coord, d_coord):
    al, ad = get_ack_ingredients(l_coord, d_coord)
    c = al
    k = 1.0 - ad
    a = 2.0 - c - k

    verify_floats([a, c, k])
    return a, c, k


def get_ack_from_ld_on_downward_face(l_coord, d_coord):
    al, ad = get_ack_ingredients(l_coord, d_coord)
    abs_k = ad
    abs_c = 1.0 - al
    abs_a = 2.0 - abs_c - abs_k

    verify_floats([abs_a, abs_c, abs_k])
    return -abs_a, -abs_c, -abs_k


def get_ack_ingredients(l_coord, d_coord, force_matrix_multiplication=False):
    if force_matrix_multiplication:
        # to reduce duplication of code
        l_q = np.array([[-l_coord], [0.0]])
        d_q = np.array([[0.0], [-d_coord]])
        Al = MATRIX_SQUARE_TO_TRIANGLE @ l_q
        Ad = MATRIX_SQUARE_TO_TRIANGLE @ d_q
        al = np.linalg.norm(Al)
        ad = np.linalg.norm(Ad)

        # because the magnitude of l_hat and d_hat is preserved as 1 by the transformation, should have no magnitude change in vectors parallel to them
        # once we know this works for sure, we can stop using matrix multiplication to get these magnitudes
        assert np.isclose(l_coord, al, atol=1e-9)
        assert np.isclose(d_coord, ad, atol=1e-9)
        return al, ad
    else:
        return float(l_coord), float(d_coord)


def get_ld_from_ack(a, c, k):
    # use `is`, not `==`, to check for positive vs negative zero, because +0.0 == -0.0
    if a is POS_ZERO:
        return (1.0, 0.0)
    elif c is POS_ZERO:
        return (0.0, 0.0)
    elif k is POS_ZERO:
        return (1.0, 1.0)
    elif a is NEG_ZERO:
        return (0.0, 1.0)
    elif c is NEG_ZERO:
        return (1.0, 1.0)
    elif k is NEG_ZERO:
        return (0.0, 0.0)
    elif is_negative_triangle_coordinate(a):
        assert is_negative_triangle_coordinate(c) and is_negative_triangle_coordinate(k), "all of a,c,k should be negative on a downward-pointing face"
        return get_ld_from_ack_on_downward_face(a, c, k)
    else:
        return get_ld_from_ack_on_upward_face(a, c, k)


def get_ld_from_ack_on_upward_face(a, c, k):
    l = c
    d = 1.0 - k
    return l, d


def get_ld_from_ack_on_downward_face(a, c, k):
    l = 1.0 - abs(c)
    d = abs(k)
    return l, d


def verify_floats(xs):
    # https://stackoverflow.com/questions/28292542/how-to-check-if-a-number-is-a-np-float64-or-np-float32-or-np-float16
    for i, x in enumerate(xs):
        assert np.issubdtype(type(x), np.floating), f"{i=}th element: {x} of type {type(x)}"


def is_negative_triangle_coordinate(x):
    return x is NEG_ZERO or x < 0


def is_positive_triangle_coordinate(x):
    return x is POS_ZERO or x > 0


def distort_ld_using_lp_transformation_in_triangle_coordinates(l, d):
    # print(f"{l=:.4f}, {d=:.4f}")
    a,c,k = get_ack_from_ld(l, d)
    # print(f"{a=:.4f}, {c=:.4f}, {k=:.4f}")
    neg = is_negative_triangle_coordinate(a)
    if neg:
        assert is_negative_triangle_coordinate(c) and is_negative_triangle_coordinate(k)
        a,c,k = -a, -c, -k
    a2 = distort.get_lp_proportion_from_theta_proportion(a)
    c2 = distort.get_lp_proportion_from_theta_proportion(c)
    k2 = distort.get_lp_proportion_from_theta_proportion(k)
    r = 2 / (a2 + c2 + k2)  # hack to try to get the new a,c,k to work with the triangular coordinates because they no longer add up to 2
    a2 *= r
    c2 *= r
    k2 *= r
    if neg:
        a2, c2, k2 = -a2, -c2, -k2
    # print(f"{a2=:.4f}, {c2=:.4f}, {k2=:.4f}")
    l2, d2 = get_ld_from_ack(a2, c2, k2)
    # print(f"{l2=:.4f}, {d2=:.4f}")
    if l2 < 0:
        assert abs(l2) < 1e-9, "negative l"
        l2 = 0.0
    if d2 < 0:
        assert abs(d2) < 1e-9, "negative d"
        d2 = 0.0
    return l2, d2



if __name__ == "__main__":
    while True:
        if random.random() < 1/4:
            # random point
            l,d = np.random.random((2,))
        elif random.random() < 1/3:
            # edge between the two faces
            r = np.random.random()
            l,d = r,r
        elif random.random() < 1/2:
            # outside edge
            a = random.choice([0.0, 1.0])
            b = random.random()
            if random.random() < 1/2:
                l,d = a,b
            else:
                l,d = b,a
        else:
            # corner
            l = random.choice([0.0, 1.0])
            d = random.choice([0.0, 1.0])
        
        l2, d2 = distort_ld_using_lp_transformation_in_triangle_coordinates(l, d)
        
        # theta_proportions = np.linspace(0, 1, 101)
        # lp_proportions = distort.get_lp_proportion_from_theta_proportion(theta_proportions)
        # diffs = lp_proportions - theta_proportions
        # assert np.isclose(diffs, -diffs[::-1], atol=1e-9).all()  # check that distortion is symmetrical on an edge (moving things closer to the center)

        print()
