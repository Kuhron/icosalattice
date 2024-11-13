# https://naturaldisasters.ai/posts/python-geopandas-world-map-tutorial/

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
# import seaborn as sns
import geopandas as gpd
import math
# import pandas as pd

import icosalattice.IcosahedronMath as icm
import icosalattice.MapCoordinateMath as mcm
import icosalattice.StartingPoints as sp


# world = gpd.read_file("https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip")
world = gpd.read_file("/home/kuhron/programming/Mapping/ne_110m_admin_0_countries.zip")

fig = plt.figure(figsize=(20, 10))
ax = fig.add_subplot()
fig.tight_layout()

# plot a basic map of the world
world.plot(
    ax=ax,
    color="lightgray",
    edgecolor="black",
    alpha=0.5
)

# add extra points and labels
starting_points, adj = sp.get_starting_points_immutable()
labels = sp.STARTING_POINT_CODES
label_to_latlon = {label: p.latlondeg() for label, p in zip(labels, starting_points)}


# double check the distances make sense so I don't make the same mistake I did before by assuming the rings were at 30 degrees of latitude
def d(a, b):
    s = 0
    for x,y in zip(a,b):
        s += (x-y)**2
    return s**0.5


for label, (lat, lon) in label_to_latlon.items():
    # https://stackoverflow.com/questions/54831344/matplotlib-plt-text-with-user-defined-circle-radii
    padding = 0.25  # in proportion of fontsize
    # ax.scatter([lon], [lat])
    ms = [-1, 1] if abs(lon) == 180 else [1]
    for m in ms:
        ax.text(m*lon, lat, label, ha="center", va="center", fontsize=20,
            bbox={"boxstyle": f"circle,pad={padding}", "fc":"yellow"})

plt.show()


