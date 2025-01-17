import numpy as np
import matplotlib.pyplot as plt


def deg_to_rad(x):
    return x * np.pi / 180


def rad_to_deg(x):
    return x * 180 / np.pi


def unit_vector_latlon_to_cartesian(lat, lon, deg=True, as_array=True):
    if deg:
        # got deg from user
        lat = deg_to_rad(lat)
        lon = deg_to_rad(lon)
    x = np.cos(lon) * np.cos(lat)
    y = np.sin(lon) * np.cos(lat)
    z = np.sin(lat)
    verify_unit_vector(x, y, z)
    if as_array:
        return np.array([x, y, z])
    else:
        return (float(x), float(y), float(z))  # get rid of np types


def unit_vector_cartesian_to_latlon(x, y, z, deg=True, as_array=True):
    # latlon [0, 0] maps to xyz [1, 0, 0] (positive x comes out of Gulf of Guinea)
    # latlon [0, 90deg] maps to xyz [0, 1, 0] (positive y comes out of Indian Ocean)
    verify_unit_vector(x, y, z)
    lat = np.arcsin(z)
    assert (abs(np.cos(lat) - np.sqrt(1 - z**2)) < 1e-6).all(), "math error in sin cos lat"
    lon = np.arctan2(y, x)  # this is the magic function I've been looking for

    if deg:
        # must give deg to user
        lat = rad_to_deg(lat)
        lon = rad_to_deg(lon)

    # input("{} {} {} -> {} {}".format(x, y, z, lat, lon))
    if as_array:
        return np.array([lat, lon])
    else:
        return (float(lat), float(lon))  # get rid of np types


def verify_unit_vector(x, y, z):
    v = np.array([x, y, z])
    assert v.shape[0] == 3
    mag = mag_3d(v)
    assert (abs(1 - mag) < 1e-6).all(), "need unit vector, but got magnitude {}\nfrom input {}".format(mag, v)


def mag_3d(v):
    assert v.shape[0] == 3  # allow underlying point array beyond this
    return np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)


def mag_3d_simple(xyz):
    x, y, z = xyz
    return (x**2 + y**2 + z**2) ** 0.5


def get_latlon_of_point_on_map(r, c, map_r_size, map_c_size,
                                map_r_min_c_min_lat, map_r_min_c_min_lon,
                                map_r_min_c_max_lat, map_r_min_c_max_lon,
                                map_r_max_c_min_lat, map_r_max_c_min_lon,
                                map_r_max_c_max_lat, map_r_max_c_max_lon,
                                deg=True):
    # print("get_latlon_of_point r={}, c={}".format(r, c))
    # like 3d printer
    # go down the rows alpha_r of the way first on left and right edge
    # then go alpha_c of the way across between those points
    alpha_r = r / map_r_size
    alpha_c = c / map_c_size
    p00 = unit_vector_latlon_to_cartesian(map_r_min_c_min_lat, map_r_min_c_min_lon, deg=deg)
    p01 = unit_vector_latlon_to_cartesian(map_r_min_c_max_lat, map_r_min_c_max_lon, deg=deg)
    p10 = unit_vector_latlon_to_cartesian(map_r_max_c_min_lat, map_r_max_c_min_lon, deg=deg)
    p11 = unit_vector_latlon_to_cartesian(map_r_max_c_max_lat, map_r_max_c_max_lon, deg=deg)

    pr0 = rotate_partially_toward_other_unit_vector(p00, p10, alpha_r)
    pr1 = rotate_partially_toward_other_unit_vector(p01, p11, alpha_r)
    prc = rotate_partially_toward_other_unit_vector(pr0, pr1, alpha_c)
    prc_latlon = unit_vector_cartesian_to_latlon(prc[0], prc[1], prc[2], deg=deg)
    return prc_latlon


def rotate_partially_toward_other_unit_vector(p, q, alpha):
    # print("rotate_partially called\n- p shape {}:\n- q shape {}:\n- alpha shape {}:".format(p.shape, q.shape, alpha.shape))
    # get vector starting from p and rotating alpha of the way towards q
    assert p.shape[0] == q.shape[0] == 3, "must apply function to 3D vectors (perhaps with larger array structure inside those three elements), got p shape {}, q shape {}".format(p.shape, q.shape)
    # assert abs(1-np.linalg.norm(p)) < 1e-6 and abs(1-np.linalg.norm(q)) < 1e-6, "must apply function to unit vectors"
    assert (0 <= alpha).all() and (alpha <= 1).all(), "alpha must be between 0 and 1"
    # if alpha == 0:
    #     return p
    # if alpha == 1:
    #     return q
    angle_p_q = angle_between_vectors(p, q)
    angle_to_move = alpha * angle_p_q
    # if angle_p_q == 0:
    #     assert p == q, "got zero angle for unequal vectors"
    #     return p
    # if angle_p_q == np.pi:
    #     assert p == -1*q, "got pi angle for non-opposite vectors"
    #     raise ValueError("great circle direction is undefined for opposite vectors")
    # https://stackoverflow.com/questions/22099490/calculate-vector-after-rotating-it-towards-another-by-angle-%CE%B8-in-3d-space
    cross = cross_3d(cross_3d(p, q), p)
    cross_mag = mag_3d(cross)
    D_tick = cross / cross_mag
    # print("angle_to_move shape {}\np shape {}\nD_tick shape {}".format(angle_to_move.shape, p.shape, D_tick.shape))
    assert p.shape[0] == 3
    assert D_tick.shape[0] == 3
    cos_array = np.cos(angle_to_move)
    sin_array = np.sin(angle_to_move)
    z_p = np.zeros((3,) + alpha.shape)
    z_d = np.zeros((3,) + alpha.shape)
    for i in range(3):
        # https://stackoverflow.com/questions/17123350/mapping-element-wise-a-numpy-array-into-an-array-of-more-dimensions
        z_p[i, ...] = cos_array * p[i]
        z_d[i, ...] = sin_array * D_tick[i]
    z = z_p + z_d
    assert z.shape[1:] == alpha.shape, "shape problem, got z shape {}".format(z.shape)
    assert z.shape[0] == 3, "shape problem, got z shape {}".format(z.shape)
    # print("- returning z")
    return z


def verify_3d_match(v1, v2):
    assert v1.shape[0] == v2.shape[0] == 3, "shape error, expected 3d vectors, perhaps with larger array structure inside those elements, got shapes {} and {}".format(v1.shape, v2.shape)
    assert v1.shape == v2.shape, "shape mismatch: {} and {}".format(v1.shape, v2.shape)


def dot_3d(v1, v2):
    verify_3d_match(v1, v2)
    point_array_shape = v1.shape[1:]
    # print("got point array shape {}".format(point_array_shape))
    res = np.zeros((1,) + point_array_shape)  # treated as single value, underlying larger array of points inside that (as though the points are many-worlds possibilities, but the array still "acts like" a single value)
    dot = np.zeros(point_array_shape)
    for i in range(3):
        dot += v1[i] * v2[i]
    return dot


def cross_3d(v1, v2):
    verify_3d_match(v1, v2)
    point_array_shape = v1.shape[1:]
    res = np.zeros((3,) + point_array_shape)  # output is a 3d vector
    a, b = v1, v2
    res[0] = a[1]*b[2] - a[2]*b[1]
    res[1] = a[2]*b[0] - a[0]*b[2]
    res[2] = a[0]*b[1] - a[1]*b[0]
    return res


def angle_between_vectors(v1, v2):
    verify_3d_match(v1, v2)
    point_array_shape = v1.shape[1:]
    dot = dot_3d(v1, v2)
    # print("dot product: {}".format(dot))
    # dot = mag(v1) * mag(v2) * cos(theta)
    len1 = mag_3d(v1)
    len2 = mag_3d(v2)
    assert (len1 > 0).all(), "v1 has zero mag: {}".format(v1)
    assert (len2 > 0).all(), "v2 has zero mag: {}".format(v2)
    cos_theta = dot / (len1 * len2)
    theta = np.arccos(cos_theta)
    # print("got theta: {}".format(theta))
    assert theta.shape == point_array_shape
    return theta


def vector_rejection_3d(v1, v2):
    # https://en.wikipedia.org/wiki/Vector_projection#Vector_rejection_2
    return v1 - (dot_3d(v1, v2) / dot_3d(v2, v2)) * v2


def xyz_distance(xyz0, xyz1):
    # make it simple and fast, no numpy
    x0, y0, z0 = xyz0
    x1, y1, z1 = xyz1
    dx = x1 - x0
    dy = y1 - y0
    dz = z1 - z0
    dxyz = (dx, dy, dz)
    return mag_3d_simple(dxyz)


def area_of_triangle_from_vertices_3d(xyz0, xyz1, xyz2):
    # https://www.quora.com/How-can-I-find-the-area-of-a-triangle-in-3D-coordinate-geometry
    x0, y0, z0 = xyz0
    x1, y1, z1 = xyz1
    x2, y2, z2 = xyz2
    ab = np.array([x1-x0, y1-y0, z1-z0])
    ac = np.array([x2-x0, y2-y0, z2-z0])
    cp = np.cross(ab, ac)
    return 1/2 * np.linalg.norm(cp)


def get_plane_containing_three_points_3d(xyz0, xyz1, xyz2):
    # plane equation of form ax*x + ay*y + az*z = c
    # https://math.stackexchange.com/questions/2686606/equation-of-a-plane-passing-through-3-points
    x0, y0, z0 = xyz0
    x1, y1, z1 = xyz1
    x2, y2, z2 = xyz2
    ab = np.array([x1-x0, y1-y0, z1-z0])
    ac = np.array([x2-x0, y2-y0, z2-z0])
    cp = np.cross(ab, ac)  # normal vector
    ax, ay, az = cp
    # now we need c, so substitute in a known point
    c = ax*x0 + ay*y0 + az*z0
    assert abs(ax*x1 + ay*y1 + az*z1 - c) < 1e-9, "error in solving for plane"
    assert abs(ax*x2 + ay*y2 + az*z2 - c) < 1e-9, "error in solving for plane"
    return ax, ay, az, c


def get_projection_dilation_ratio_of_point_onto_plane(xyz, ax, ay, az, c):
    # plane equation of form ax*x + ay*y + az*z = c
    xp, yp, zp = xyz
    c_of_p = ax*xp + ay*yp + az*zp
    if abs(c_of_p) < 1e-9:
        ratio = 1 if abs(c) < 1e-9 else float("inf") if c > 0 else -float("inf")
    else:
        ratio = c/c_of_p
    return ratio


def project_point_onto_plane(xyz, ax, ay, az, c):
    # plane equation of form ax*x + ay*y + az*z = c
    r = get_projection_dilation_ratio_of_point_onto_plane(xyz, ax, ay, az, c)
    return xyz * r


def get_unit_sphere_midpoint_from_latlon(latlon0, latlon1, as_array=True):
    xyz0 = unit_vector_latlon_to_cartesian(*latlon0)
    xyz1 = unit_vector_latlon_to_cartesian(*latlon1)
    xyz = get_unit_sphere_midpoint_from_xyz(xyz0, xyz1)
    return unit_vector_cartesian_to_latlon(*xyz, as_array=as_array)


def get_unit_sphere_midpoint_from_xyz(xyz0, xyz1, as_array=True):
    # do it simple with basic math functions, no numpy casting or anything fancy, want fast
    x0, y0, z0 = xyz0
    x1, y1, z1 = xyz1
    xm = (x0 + x1) / 2
    ym = (y0 + y1) / 2
    zm = (z0 + z1) / 2
    m_raw = (xm, ym, zm)
    mag = mag_3d_simple(m_raw)
    m = (xm / mag, ym / mag, zm / mag)
    assert abs(mag_3d_simple(m) - 1) < 1e-9
    if as_array:
        return np.array(m)
    else:
        return m


def get_xyz_from_latlon(latlon, deg=True, as_array=True):
    # alias function
    return unit_vector_latlon_to_cartesian(*latlon, deg=deg, as_array=as_array)


def get_latlon_from_xyz(xyz, deg=True, as_array=True):
    # alias function
    return unit_vector_cartesian_to_latlon(*xyz, deg=deg, as_array=as_array)



# ---- UNSORTED STUFF BELOW ---- I have quarantined it so that I can draw on it if needed in icosalattice library, but won't be splitting files into part that lives in this repo and the rest that I don't need living in another file outside the repo ---- #

"""

def get_radius_about_center_surface_point_for_circle_of_area_proportion_on_unit_sphere(area_sphere_proportion):
    assert 0 < area_sphere_proportion < 1, "expected size must be proportion of sphere surface area between 0 and 1, but got {}".format(area_sphere_proportion)
    # expected size is in terms of proportion of sphere surface area; note that this does not scale linearly with radius in general
    # because will be using Euclidean distance in R^3, need to do some trig to convert the surface area proportion to 3d radius
    # center point is on the unit sphere, if r=1, what is the area within that distance? (the distance is a chord through the sphere's interior), the total sphere surface area = 4*pi*r^2 = 4*pi
    # f(r) = integral(0 to r, dA/dr dr); f(2) = the whole sphere = 4*pi; f(sqrt(2)) = half sphere = 2*pi
    # dA = 2*pi*r' dr, where r' is the radius of the flat circle that r points to, from the central axis which runs through the starting point and the sphere's center
    # drew pictures and got that r'^2 = r^2 - r^4/4; checked r'(r=sqrt(2))=1, r'(r=0)=0, r'(r=2)=0, r'(r=1)=sqrt(3)/2, all work
    # so dA = 2*pi*sqrt(r^2 - r^4/4) dr
    # but dA should be scaled up by some trig factor (e.g. it will be sqrt(2) times greater when it is slanted at 45 deg, and 0 times greater when it is vertical)
    # dA/dr' = 2*pi*dl, imagine lowering the circle by dh, so that r' rises by dr', then the slanted line on the sphere's surface is dl
    # if theta is angle with vertical axis, cos theta = dr'/dl, so dl = dr'/cos(theta), and from the center, see that sin(theta) r'/1
    # so cos(theta) = sqrt(1-r'^2), so dA = 2*pi* r' * dr'/sqrt(1-r'^2)
    # can integrate over r' instead of r now
    # f(r) = 2*pi* integral(0 to sqrt(r^2 - r^4/4), r'/sqrt(1-r'^2) dr')
    # proportion(r) = 1/(4*pi) * f(r)
    # integral(s/sqrt(1-s^2) ds) = -1*sqrt(1-s^2) => (lots of whiteboard scribbles) => proportion(r) = r^2/4 (for r in [0, 2])
    # => r(proportion) = 2*sqrt(proportion)
    radius_from_center_in_3d = 2 * np.sqrt(area_sphere_proportion)
    return radius_from_center_in_3d

    
"""