# trying to figure out how to fix the ack sum to be 2 on both adjust and readjust of ld so that the inverse transformation works correctly

import numpy as np
import matplotlib.pyplot as plt

from icosalattice.TriangularPeelCoordinates import get_ack_from_ld
import icosalattice.FacePlaneDistortion as distort

vals = np.linspace(0, 1, 500)
lds = [(l,d) for d in vals for l in vals]
sums = []
for l,d in lds:
    a,c,k = get_ack_from_ld(l, d)
    neg = a<0 or c<0 or k<0
    if neg:
        assert a<=0 and c<=0 and k<=0, (a,c,k)
        a,c,k = -a,-c,-k
    
    a2 = distort.get_lp_proportion_from_theta_proportion(a)
    c2 = distort.get_lp_proportion_from_theta_proportion(c)
    k2 = distort.get_lp_proportion_from_theta_proportion(k)

    if neg:
        a2,c2,k2 = -a2,-c2,-k2

    s = abs(a2 + c2 + k2)
    if abs(s-2) > 0.1:
        print(f"{l=:.6f}, {d=:.6f}\n{a=:.6f}, {c=:.6f}, {k=:.6f}\n{a2=:.6f}, {c2=:.6f}, {k2=:.6f}")
        input("check")
    sums.append(s)

print(min(sums), max(sums))
plt.hist(sums, bins=100)
plt.show()
