# https://naturaldisasters.ai/posts/python-geopandas-world-map-tutorial/

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
# import seaborn as sns
import geopandas as gpd
import math
# import pandas as pd

import icosalattice.IcosahedronMath as icm
import icosalattice.MapCoordinateMath as mcm


world = gpd.read_file("https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip")

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

# turn off axis ticks
# ax.set_xticks([])
# ax.set_yticks([])

# plt.title("Basic Map of World with GeoPandas")

# add extra points and labels
label_to_latlon = {
    "A": (90, 0),
    "B": (-90, 0),
    "C": (icm.RING_LAT_DEG, 0),
    "D": (-icm.RING_LAT_DEG, 36),
    "E": (icm.RING_LAT_DEG, 72),
    "F": (-icm.RING_LAT_DEG, 108),
    "G": (icm.RING_LAT_DEG, 144),
    "H": (-icm.RING_LAT_DEG, 180),
    "I": (icm.RING_LAT_DEG, -144),
    "J": (-icm.RING_LAT_DEG, -108),
    "K": (icm.RING_LAT_DEG, -72),
    "L": (-icm.RING_LAT_DEG, -36),
}

# double check the distances make sense so I don't make the same mistake I did before by assuming the rings were at 30 degrees of latitude
def d(a, b):
    s = 0
    for x,y in zip(a,b):
        s += (x-y)**2
    return s**0.5

# letters = "ABCDEFGHIJKL"
# ds = {letters[i]: {letters[j]: None for j in range(12)} for i in range(12)}
# for i in range(12):
#     for j in range(12):
#         lli = label_to_latlon[letters[i]]
#         llj = label_to_latlon[letters[j]]
#         xyzi = mcm.unit_vector_latlon_to_cartesian(lli[0], lli[1], deg=True)
#         xyzj = mcm.unit_vector_latlon_to_cartesian(llj[0], llj[1], deg=True)
#         ds[letters[i]][letters[j]] = d(xyzi, xyzj)
# df = (pd.DataFrame.from_dict(ds)*1000).round().astype(int)
# print(df)

for label, (lat, lon) in label_to_latlon.items():
    # https://stackoverflow.com/questions/54831344/matplotlib-plt-text-with-user-defined-circle-radii
    padding = 0.25  # in proportion of fontsize
    # ax.scatter([lon], [lat])
    ms = [-1, 1] if abs(lon) == 180 else [1]
    for m in ms:
        ax.text(m*lon, lat, label, ha="center", va="center", fontsize=20,
            bbox={"boxstyle": f"circle,pad={padding}", "fc":"yellow"})

plt.show()


