# trying to figure out how to fix the ack sum to be 2 on both adjust and readjust of ld so that the inverse transformation works correctly
# if can get a function for what the sum factor is at each point then we can just divide by it and then do the inversion of lp
# expect it to obey the triangle's symmetry

import numpy as np
import matplotlib.pyplot as plt

from icosalattice.TriangularPeelCoordinates import get_ack_from_ld
import icosalattice.FacePlaneDistortion as distort

ALPHA = distort.ALPHA

n = 100
vals = np.linspace(0, 1, n)
X = np.zeros((n,n))
Y = np.zeros((n,n))
sums = np.zeros((n, n))

lp = distort.get_lp_proportion_from_theta_proportion
lp_inv = distort.get_theta_proportion_from_lp_proportion

MIN_LP_SUM = 3 * (1/2 + (1/4 + np.sqrt(5)/4) * np.tan(ALPHA/6))
MAX_LP_INV_SUM = 3 * (1/2 + np.atan(1/6 * (np.sqrt(5) - 1)) / ALPHA)

R32 = 3**0.5 / 2
xA, yA = 0, R32 * distort.W
xC, yC = 0.5 * distort.W, 0
xL, yL = 0, -R32 * distort.W


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

        # a2 = lp(a)
        # c2 = lp(c)
        # k2 = lp(k)
        a2 = lp_inv(a)
        c2 = lp_inv(c)
        k2 = lp_inv(k)

        if neg:
            a2,c2,k2 = -a2,-c2,-k2

        s = abs(a2 + c2 + k2)
        if abs(s-2) > 0.1:
            print(f"{l=:.6f}, {d=:.6f}\n{a=:.6f}, {c=:.6f}, {k=:.6f}\n{a2=:.6f}, {c2=:.6f}, {k2=:.6f}")
            input("check")
        sums[i,j] = s

assert np.isclose(3 * lp(2/3), MIN_LP_SUM, rtol=1e-9)
assert np.isclose(3 * lp_inv(2/3), MAX_LP_INV_SUM, rtol=1e-9)

# print(MIN_LP_SUM * MAX_LP_INV_SUM / 4)  # they are not reciprocals

plt.pcolormesh(X, Y, sums)
plt.colorbar()
plt.gca().set_aspect("equal")
plt.show()

