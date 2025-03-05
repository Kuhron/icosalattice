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

# discoveries while trying to figure out how to get r1 = 2/lp_sum_orig as a function of (a2, c2, k2) (the post-adjustment point)
# # because I do think this is a function and this would be the clearest way to invert the adjustment
# - r1 is not a function of geometric mean of (a2, c2, k2) (multiple r1 values for one geomean value)
# - 

# alternate idea similar to cpg1:
# - draw shortest lines from the point to each side, lp those intersection points, see where the perpendiculars from the adjusted edge points get you
# - very likely they will not converge, so like with cpg1 you'll have to adjust them somehow
# - maybe the point that minimizes sum of distances to those three?


import numpy as np
import matplotlib.pyplot as plt

from icosalattice.TriangularPeelCoordinates import get_ack_from_ld, get_ld_from_ack
import icosalattice.FacePlaneDistortion as distort
from icosalattice.MathUtil import zigzag, mod


n = 500
vals = sorted(list(np.linspace(0, 1, n)) + [1/3, 2/3])
n = len(vals)


ALPHA = distort.ALPHA

MIN_LP_SUM = 3 * (1/2 + (1/4 + np.sqrt(5)/4) * np.tan(ALPHA/6))
MAX_LP_INV_SUM = 3 * (1/2 + np.atan(1/6 * (np.sqrt(5) - 1)) / ALPHA)

R32 = 3**0.5 / 2
xA, yA = 0, R32 * distort.W
xC, yC = 0.5 * distort.W, 0
xK, yK = -xC, yC
xL, yL = 0, -R32 * distort.W
xQ, yQ = 1/3*(xA+xC+xK), 1/3*(yA+yC+yK)  # the center of the upward-pointing face (ACK)

# original x,y of points
X1 = np.zeros((n,n))
Y1 = np.zeros((n,n))

# x,y of points after adjustment
X2 = np.zeros((n,n))
Y2 = np.zeros((n,n))

lp_sums = np.zeros((n, n))
lp_inv_sums = np.zeros((n, n))

rs = np.zeros((n, n))  # what value of r = 2/lp_sum_orig got us to this point

rhos_orig = np.zeros((n,n))  # radius from center of face, before adjusting the point
rhos_adj = np.zeros((n,n))  # radius from center of face, after adjusting the point
rho_rels_orig = np.zeros((n,n))  # radius as proportion of max radius at that angle, before adjusting the point
rho_rels_adj = np.zeros((n,n))  # radius as proportion of max radius at that angle, after adjusting the point

thetas_orig = np.zeros((n,n))  # angle on one-sixth triangle face (WLOG so that center of edge is at theta=0 and corner of triangle is at theta=pi/3), before adjusting the point
thetas_adj = np.zeros((n,n))  # theta after adjusting the point
theta_rels_orig = np.zeros((n,n))  # theta as proportion of max theta (=pi/3), before adjusting the point
theta_rels_adj = np.zeros((n,n))  # theta proportion after adjusting the point

lp = distort.get_lp_proportion_from_theta_proportion
lp_inv = distort.get_theta_proportion_from_lp_proportion

get_rho_from_xy = lambda x,y: ((x-xQ)**2 + (y-yQ)**2) ** 0.5
get_theta_from_xy = lambda x,y: theta_zig(get_theta_0_from_xy(x,y))
get_theta_0_from_xy = lambda x,y: theta_shift(get_theta_raw_from_xy(x,y))
get_theta_raw_from_xy = lambda x,y: 0 if xy_is_centroid(x,y) else np.arctan2((y-yQ),(x-xQ))
theta_shift = lambda theta: theta - np.pi/6
theta_zig = lambda theta: zigzag(mod(theta, 2*np.pi), np.pi/3)

rho_max = 2/3 * distort.B
rho_min = 1/3 * distort.B
get_rho_max_for_theta = lambda theta: rho_min / np.cos(theta)  # maximum distance (or radius) from centroid on face plane
get_rho_relative = lambda rho, theta: rho / get_rho_max_for_theta(theta)
get_theta_relative = lambda theta: theta / (np.pi/3)
xy_is_centroid = lambda x,y: np.isclose(x, xQ, rtol=1e-12) and np.isclose(y, yQ, rtol=1e-12)


for i in range(n):
    for j in range(n):
        l1 = vals[i]
        d1 = vals[j]

        x1 = xC + l1 * (xA - xC) + d1 * (xL - xC)
        y1 = yC + l1 * (yA - yC) + d1 * (yL - yC)
        X1[i,j] = x1
        Y1[i,j] = y1

        a1,c1,k1 = get_ack_from_ld(l1, d1)
        neg = a1<0 or c1<0 or k1<0
        if neg:
            assert a1<=0 and c1<=0 and k1<=0, (a1,c1,k1)
            a1,c1,k1 = -a1,-c1,-k1

        a2_orig = lp(a1)
        c2_orig = lp(c1)
        k2_orig = lp(k1)
        lp_sum = a2_orig + c2_orig + k2_orig
        r1 = 2 / lp_sum
        a2 = r1 * a2_orig
        c2 = r1 * c2_orig
        k2 = r1 * k2_orig
        # crucially, the inverse sum that we can access during deadjustment is the sum of inverses at the values that we got *after* adjustment!
        a3_orig = lp_inv(a2)
        c3_orig = lp_inv(c2)
        k3_orig = lp_inv(k2)
        lp_inv_sum = a3_orig + c3_orig + k3_orig
        r2 = 2 / lp_inv_sum
        a3 = r2 * a3_orig
        c3 = r2 * c3_orig
        k3 = r2 * k3_orig

        # # for printing to Excel-readable file
        print(f"{a1:.6f} {c1:.6f} {k1:.6f} {a1+c1+k1:.6f} {a2_orig:.6f} {c2_orig:.6f} {k2_orig:.6f} {lp_sum:.6f} {a2:.6f} {c2:.6f} {k2:.6f} {a2+c2+k2:.6f}")

        # for printing more human-readably for me to think about what to do
        # print(f"({a:.6f}, {c:.6f}, {k:.6f}) (sum {a+c+k:.6f}) >lp> ({a2_orig:.6f}, {c2_orig:.6f}, {k2_orig:.6f}) (sum {lp_sum_orig:.6f}) >r> ({a2:.6f}, {c2:.6f}, {k2:.6f}) (sum {a2+c2+k2:.6f})")

        if neg:
            a2,c2,k2 = -a2,-c2,-k2

        l2, d2 = get_ld_from_ack(a2, c2, k2)
        x2 = xC + l2 * (xA - xC) + d2 * (xL - xC)
        y2 = yC + l2 * (yA - yC) + d2 * (yL - yC)

        X2[i,j] = x2
        Y2[i,j] = y2

        if lp_sum < MIN_LP_SUM - 1e-9 or lp_sum > 2 + 1e-9 or lp_inv_sum > MAX_LP_INV_SUM + 1e-9 or lp_inv_sum < 2 - 1e-9:
            raise Exception(f"{l1=:.6f}, {d1=:.6f}\n{a1=:.6f}, {c1=:.6f}, {k1=:.6f}\n{a2=:.6f}, {c2=:.6f}, {k2=:.6f}")
        lp_sums[i,j] = lp_sum
        lp_inv_sums[i,j] = lp_inv_sum
        rs[i, j] = r1

        rho_orig = get_rho_from_xy(x1, abs(y1))
        rho_adj = get_rho_from_xy(x2, abs(y2))
        theta_orig = get_theta_from_xy(x1, abs(y1))
        theta_adj = get_theta_from_xy(x2, abs(y2))

        rhos_orig[i, j] = rho_orig
        rhos_adj[i, j] = rho_adj
        rho_rels_orig[i, j] = get_rho_relative(rho_orig, theta_orig)
        rho_rels_adj[i, j] = get_rho_relative(rho_adj, theta_adj)
        thetas_orig[i, j] = theta_orig
        thetas_adj[i, j] = theta_adj
        theta_rels_orig[i, j] = get_theta_relative(theta_orig)
        theta_rels_adj[i, j] = get_theta_relative(theta_adj)


assert np.isclose(3 * lp(2/3), MIN_LP_SUM, rtol=1e-9)
assert np.isclose(3 * lp_inv(2/3), MAX_LP_INV_SUM, rtol=1e-9)

# print(f"{rho_max = }, {rhos_orig.max() = }")
# # if rho_max is ~0.577, you are using the wrong triangle size (side length 1)
# # you should see rho_max = ~0.607, when side length = W ~= 1.05


# for arr in [thetas_orig, thetas_adj]:
#     plt.pcolormesh(X2, Y2, arr, cmap="jet")
#     plt.colorbar()
#     plt.gca().set_aspect("equal")
#     plt.show()

plt.scatter(rho_rels_orig, theta_rels_orig, c=rho_rels_adj-rho_rels_orig, cmap="jet")
plt.xlabel("relative rho of original point")
plt.ylabel("relative theta of original point")
plt.colorbar().set_label("difference in relative rho: adjusted point minus original point")
plt.show()

plt.scatter(rho_rels_orig, theta_rels_orig, c=theta_rels_adj-theta_rels_orig, cmap="jet")
plt.xlabel("relative rho of original point")
plt.ylabel("relative theta of original point")
plt.colorbar().set_label("difference in relative theta: adjusted point minus original point")
plt.show()
