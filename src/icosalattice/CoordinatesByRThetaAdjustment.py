import math
import numpy as np

import icosalattice.PeelCoordinates as pe
from icosalattice.FacePlaneDistortion import ALPHA, B, G, GAMMA, H, W
from icosalattice.MathUtil import zigzag, zigzag_inverse, mod, round_off_unwanted_float_precision
import icosalattice.StartingPoints as sp
from icosalattice.CoordinatesByAncestry import get_xyz_of_initial_point_code


R3 = 3**0.5
R5 = 5**0.5

R = 1  # sphere radius
# GAMMA is the angle from perspective of sphere center from centroid of face plane to vertex

x_A = 0
y_A = 2/3 * R3/2 * W

x_C = 1/2 * W
y_C = -1/3 * R3/2 * W

x_K = -x_C
y_K = y_C

x_L = 0
y_L = -4/3 * R3/2 * W

pA = (x_A, y_A)
pC = (x_C, y_C)
pK = (x_K, y_K)
pL = (x_L, y_L)
pQ = (0, 0)
pAC = (1/2 * (x_A + x_C), 1/2 * (y_A + y_C))

DXL = x_A - x_C  # dx due to one unit of l-ward motion (L = "left" in peel coordinates, not the point L)
DYL = y_A - y_C
DXD = x_L - x_C  # dx due to one unit of d-ward motion (D = "down" in peel coordinates, not the point D)
DYD = y_L - y_C

# mapping from (l,d) to (x,y) offset from pC
M = np.array([
    [DXL, DXD],
    [DYL, DYD],
])
M_inv = np.linalg.inv(M)
(DLX, DLY), (DDX, DDY) = M_inv



def get_xyz_from_point_code_using_r_theta_adjustment(pc, as_array=True):
    if pc in sp.STARTING_POINT_CODES:
        return get_xyz_of_initial_point_code(pc)
    
    spc, l, d = pe.get_raw_peel_coordinates_from_point_code(pc)

    # put all points on an upward-facing triangle; the (l,d) transformation respects all symmetries of an equilateral triangle
    # so we can just flip it over the line in the 2 direction from the ancestral point
    if d > l:
        flipped = True
        l_wlog = d
        d_wlog = l
    else:
        flipped = False
        l_wlog = l
        d_wlog = d
    
    l_new_wlog, d_new_wlog = transform_ld_by_r_theta_adjustment(l_wlog, d_wlog)

    if flipped:
        l_new = d_new_wlog
        d_new = l_new_wlog
    else:
        l_new = l_new_wlog
        d_new = d_new_wlog
    
    return pe.get_xyz_from_adjusted_peel_coordinates((spc, l_new, d_new), as_array=as_array)


def transform_ld_by_r_theta_adjustment(l, d):
    # copying a bunch of stuff from my Desmos graph "TriangleVectorField": https://www.desmos.com/calculator/qydlzra4qu

    # map all faces onto an upward-pointing equilateral triangle centered on the origin, with edge length w
    # within that, map all points into one-sixth of the triangle (between the centroid, midpoint of CA, and A)

    # transformations from interval [0,1] to itself
    # lp = lambda x: (np.sin(ALPHA/2) - np.cos(ALPHA/2)*np.tan(ALPHA/2 - x*ALPHA)) / (2*np.sin(ALPHA/2))
    lq = lambda x: 1/2 * (1 + R5) * np.tan(1/2 * x * ALPHA)

    f_theta_zig = lambda theta: zigzag(mod(theta, 2*np.pi), a=np.pi/3)
    f_theta_unzig = lambda theta, theta_0: zigzag_inverse(theta, a=np.pi/3, n=math.floor(theta_0 / (np.pi/3)))

    f_rho = lambda x, y: (x**2 + y**2)**0.5  # distance (or radius) from centroid on face plane
    rho_max = 2/3 * B
    rho_min = 1/3 * B
    f_rho_max_for_theta = lambda theta: rho_min / np.cos(theta)  # maximum distance (or radius) from centroid on face plane
    # f_rho_relative = lambda rho, theta: rho / f_rho_max_for_theta(theta)  # proportion of maximum possible distance from centroid on face plane
    f_sigma_new_from_centroid = lambda rho: GAMMA * rho / rho_max  # angle that is the same proportion of gamma as the rho is of the rho to a vertex on the face plane
    # f_s_for_rho_raw = lambda rho: R * f_sigma_new_from_centroid(rho)  # arc length from centroid on sphere surface such that this arc length is the same proportion of the arc length to a vertex as the rho is of the rho to a vertex on the face plane
    f_rho_2 = lambda rho: G * np.tan(f_sigma_new_from_centroid(rho))

    # get the adjusted new radius by scaling linearly so that the maximum radius at old theta maps to the maximum radius at new theta
    f_rho_3 = lambda rho, theta_old, theta_new: f_rho_2(rho) * f_rho_max_for_theta(theta_new) / f_rho_2(f_rho_max_for_theta(theta_old))

    f_theta_shift = lambda theta: theta - np.pi/6  # set angle of AC midpoint to 0 and angle of A to 60 deg
    f_theta_unshift = lambda theta: theta + np.pi/6
    # theta_max = 60 * np.pi/180  # maximum angle from perspective of centroid from 0-degree reference to the point (found at the vertices)

    # get new adjusted angle for the point
    f_theta_2 = lambda theta: np.atan(f_j_2(theta) / rho_min)

    j_max = W/2  # maximum length along AC edge for any point (found at the vertices e.g. A)
    f_j_theta = lambda theta: rho_min * np.tan(theta)  # j (displacement from AC midpoint along AC edge) for a given theta away from AC from perspective of centroid
    f_j_relative = lambda theta: f_j_theta(theta) / j_max
    f_j_2 = lambda theta: j_max * lq(f_j_relative(theta))
    
    # now to actually use the given point
    # print(f"\n{l = :f}, {d = :f}")
    dx = l * DXL + d * DXD
    dy = l * DYL + d * DYD
    x = x_C + dx
    y = y_C + dy

    rho = f_rho(x, y)
    theta_raw = np.atan2(y, x)
    theta_from_AC = f_theta_shift(theta_raw)
    theta = f_theta_zig(theta_from_AC)
    theta_2 = f_theta_2(theta)
    # rho_relative = f_rho_relative(rho, theta)
    # rho_max_for_theta = f_rho_max_for_theta(theta)
    # j_theta = f_j_theta(theta)
    # j_relative = f_j_relative(theta)
    # j_2 = f_j_2(theta)
    theta_new = f_theta_unzig(theta_2, theta_from_AC)
    # rho_2 = f_rho_2(rho)
    rho_new = f_rho_3(rho, theta, theta_2)
    
    # print(f"{rho = :f} -> {rho_new:f}")
    # print(f"theta (from AC) = {theta_from_AC/np.pi:f} pi -> {theta_new/np.pi:f} pi")
    assert -1e-9 <= rho_new <= rho_max + 1e-9, f"expected 0 <= {rho_new = :f} <= {rho_max:f}"

    x_new = rho_new * np.cos(f_theta_unshift(theta_new))
    y_new = rho_new * np.sin(f_theta_unshift(theta_new))
    dx_new = x_new - x_C
    dy_new = y_new - y_C
    # print(f"{x = :f} -> {x_new:f}\n{y = :f} -> {y_new:f}")

    l_C, d_C = 0, 0
    dl_new = dx_new * DLX + dy_new * DLY
    dd_new = dx_new * DDX + dy_new * DDY
    l_new = l_C + dl_new
    d_new = d_C + dd_new

    l_new = round_off_unwanted_float_precision(l_new)
    d_new = round_off_unwanted_float_precision(d_new)

    # print(f"{l = :f} -> {l_new:f}\n{d = :f} -> {d_new:f}\n")
    return l_new, d_new



if __name__ == "__main__":
    for pc in [
        "C", "C1", "C2", "C11", "C01", "C22", "C02",
        "A", "B", "D", "K",
    ]:
        xyz = get_xyz_from_point_code_using_r_theta_adjustment(pc, as_array=False)
        print(xyz)
        input("check")
