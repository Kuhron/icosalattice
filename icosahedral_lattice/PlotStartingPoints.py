import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd

from icosahedral_lattice.IcosahedronMath import STARTING_POINTS

def plot_starting_points():

    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot()

    # plot a basic map of the world
    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    world.plot(
        ax=ax,
        color="lightgray",
        edgecolor="black",
        alpha=0.5
    )
