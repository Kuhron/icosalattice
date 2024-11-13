import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely import transform
import math


# TODO figure out how to split up the world into the five peels and plot the portion of the map found on each triangular face


def sign(x):
    return 1 if x > 0 else -1 if x < 0 else 0

def transform_coords(lonlat):
    assert type(lonlat) is np.ndarray
    N, two = lonlat.shape
    assert two == 2, lonlat.shape
    return np.apply_along_axis(transform_one_coord_pair, 1, lonlat)

def transform_one_coord_pair(lonlat):
    lon, lat = lonlat
    new_lon = -lon + np.sin(lat)
    new_lat = lat + lon**2 * 70/(180**2)
    return (new_lon, new_lat)


if __name__ == "__main__":
    # world = gpd.read_file("https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip")
    world = gpd.read_file("/home/kuhron/programming/Mapping/ne_110m_admin_0_countries.zip")
    for i in world.index:
        world.loc[i, "geometry"] = transform(world.loc[i, "geometry"], transformation=transform_coords)

    world.plot(
        ax=plt.gca(),
        color="lightgray",
        edgecolor="black",
        alpha=0.5
    )
    plt.gcf().tight_layout()
    plt.gca().set_aspect("equal")
    plt.show()

