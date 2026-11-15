import icosalattice.MapCoordinateMath as mcm



def get_point_code_from_xyz(xyz, method, as_array=True):
    raise NotImplementedError("do some kind of binary search thing")


def get_point_code_from_latlon(latlon, method, as_array=True):
    xyz = mcm.get_xyz_from_latlon(latlon, as_array=as_array)
    return get_point_code_from_xyz(xyz, method, as_array=as_array)


if __name__ == "__main__":
    latlons = [
        [47.37444, 8.54111],  # Zurich, Switzerland
        [10.77556, 106.70194],  # Ho Chi Minh City, Vietnam
        [-54.807222, -68.308056],  # Ushuaia, Argentina
        [44.052222, -123.0925],  # Eugene, Oregon, United States
    ]

    raise NotImplementedError