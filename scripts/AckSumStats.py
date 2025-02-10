# trying to figure out how to fix the ack sum to be 2 on both adjust and readjust of ld so that the inverse transformation works correctly
# if can get a function for what the sum factor is at each point then we can just divide by it and then do the inversion of lp
# expect it to obey the triangle's symmetry

# what I really need is a one-to-one correspondence between the lp sum and the lp inverse sum
# - so that after I've lp-inversed the points for deadjustment and summed them, I know which lp factor was used to get here so I can divide by it
# don't actually need lp sum or lp inverse sum as function of point, just need correspondence between these two value series and require one-to-one

# discovery: lp_sum_orig and lp_inv_sum_orig are NOT in one-to-one correspondence!
# - hopefully there is some other quantity we can get from the adjusted ack to know what lp_sum was used to adjust them
# - maybe the hack of multiplying by 2/lp_sum_orig is not one-to-one on the triangle? should check this, hope it's not folding the fabric over itself
# - - from Desmos it looks okay by moving the test point around, but beware of this possibility just in case

# discovery from comparing Desmos screenshots: the contour lines on the lp_sums and lp_inv_sums are NOT the same!
# - they all look like rounded triangles, but the counterpart contour does not have the same shape
# - that is, if you take the lp_sum contour passing through a given point on the triangle, and compare the lp_inv_sum contour through the same point,
# - they will both be rounded triangles but one of them will bulge out more than the other
# - to me this further shows that multiplying by the factor 2/lp_sum_orig is too simple and leads to difficulty in inverting the adjustment

# alternate idea similar to cpg1:
# - draw shortest lines from the point to each side, lp those intersection points, see where the perpendiculars from the adjusted edge points get you
# - very likely they will not converge, so like with cpg1 you'll have to adjust them somehow
# - maybe the point that minimizes sum of distances to those three?


import numpy as np
import matplotlib.pyplot as plt

from icosalattice.TriangularPeelCoordinates import get_ack_from_ld
import icosalattice.FacePlaneDistortion as distort

ALPHA = distort.ALPHA

n = 500
vals = sorted(list(np.linspace(0, 1, n)) + [1/3, 2/3])
n = len(vals)

X = np.zeros((n,n))
Y = np.zeros((n,n))
sums = np.zeros((n, n))

rs_at_adjusted_points = np.zeros((n, n))  # what value of r = 2/lp_sum_orig got us to this point

lp = distort.get_lp_proportion_from_theta_proportion
lp_inv = distort.get_theta_proportion_from_lp_proportion

MIN_LP_SUM = 3 * (1/2 + (1/4 + np.sqrt(5)/4) * np.tan(ALPHA/6))
MAX_LP_INV_SUM = 3 * (1/2 + np.atan(1/6 * (np.sqrt(5) - 1)) / ALPHA)

R32 = 3**0.5 / 2
xA, yA = 0, R32 * distort.W
xC, yC = 0.5 * distort.W, 0
xL, yL = 0, -R32 * distort.W

lp_sums = []
lp_inv_sums = []

for i in range(n):
    for j in range(n):
        l = vals[i]
        d = vals[j]

        x = xC + l * (xA - xC) + d * (xL - xC)
        y = yC + l * (yA - yC) + d * (yL - yC)
        X[i,j] = x
        Y[i,j] = y

        a,c,k = get_ack_from_ld(l, d)
        neg = a<0 or c<0 or k<0
        if neg:
            assert a<=0 and c<=0 and k<=0, (a,c,k)
            a,c,k = -a,-c,-k

        a2_orig = lp(a)
        c2_orig = lp(c)
        k2_orig = lp(k)
        lp_sum_orig = a2_orig + c2_orig + k2_orig
        r1 = 2 / lp_sum_orig
        a2 = r1 * a2_orig
        c2 = r1 * c2_orig
        k2 = r1 * k2_orig
        # crucially, the inverse sum that we can access during deadjustment is the sum of inverses at the values that we got *after* adjustment!
        a3_orig = lp_inv(a2)
        c3_orig = lp_inv(c2)
        k3_orig = lp_inv(k2)
        lp_inv_sum_orig = a3_orig + c3_orig + k3_orig
        r2 = 2 / lp_inv_sum_orig
        a3 = r2 * a3_orig
        c3 = r2 * c3_orig
        k3 = r2 * k3_orig

        lp_sums.append(lp_sum_orig)
        lp_inv_sums.append(lp_inv_sum_orig)

        # # for printing to Excel-readable file
        # print(f"{a:.6f} {c:.6f} {k:.6f} {a+c+k:.6f} {a2_orig:.6f} {c2_orig:.6f} {k2_orig:.6f} {lp_sum_orig:.6f} {a3_orig:.6f} {c3_orig:.6f} {k3_orig:.6f} {lp_inv_sum_orig:.6f}")

        # for printing more human-readably for me to think about what to do
        print(f"({a:.6f}, {c:.6f}, {k:.6f}) (sum {a+c+k:.6f}) >lp> ({a2_orig:.6f}, {c2_orig:.6f}, {k2_orig:.6f}) (sum {lp_sum_orig:.6f}) >r> ({a2:.6f}, {c2:.6f}, {k2:.6f}) (sum {a2+c2+k2:.6f})")

        if neg:
            a2,c2,k2 = -a2,-c2,-k2

        s = lp_sum_orig
        # s = lp_inv_sum_orig

        if s < MIN_LP_SUM - 1e-9 or s > MAX_LP_INV_SUM + 1e-9:
            raise Exception(f"{l=:.6f}, {d=:.6f}\n{a=:.6f}, {c=:.6f}, {k=:.6f}\n{a2=:.6f}, {c2=:.6f}, {k2=:.6f}")
        sums[i,j] = s
        rs_at_adjusted_points[i, j] = r1

assert np.isclose(3 * lp(2/3), MIN_LP_SUM, rtol=1e-9)
assert np.isclose(3 * lp_inv(2/3), MAX_LP_INV_SUM, rtol=1e-9)

# print(MIN_LP_SUM * MAX_LP_INV_SUM / 4)  # they are not reciprocals

# plt.pcolormesh(X, Y, sums, cmap="jet")
plt.pcolormesh(X, Y, rs_at_adjusted_points, cmap="jet")
plt.colorbar()
plt.gca().set_aspect("equal")
plt.show()

# plt.scatter(lp_sums, lp_inv_sums)
# plt.show()
