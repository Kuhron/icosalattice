# icosalattice
This is a Python library which implements a lattice on the sphere (by this I mean the 2-sphere, which is the surface of a 3D sphere). The lattice is based on subdividing the edges and faces of an icosahedron to arbitrary precision.

## What does it do? / What is it for?
The lattice provides a coordinate system for specifying points on the sphere. It is different from either 3D Cartesian coordinates or spherical coordinates. 

Each point has a code, which is a string telling exactly where it is based on how it was generated from points in the previous iteration. 

The iterative way in which the lattice is constructed also creates a natural graph of connections between points. There is a complicated but efficient "arithmetic" for determining the neighbors of a point based on its code. This makes it useful for simulating things spreading around on the sphere.

## Why do this?
This lattice was originally conceived as a way to generate terrain for a made-up planet. I wanted a set of points roughly evenly spaced across the sphere. However, due to Gauss's Theorema Egregium, we cannot use lattices in latitude-longitude space for this because they will oversample near the poles. I decided to take the Platonic solid which is closest to a sphere and use that as the basis for a coordinate system. Since the faces of the icosahedron are triangles, this also results in all points except the original 12 vertices having 6 neighbors. Any area of the globe not including one of these vertices resembles a hexagonal/triangular lattice, which I find nicer than a square lattice.

# How the lattice works

Place an icosahedron just inside the unit sphere such that one of its vertices lies at the sphere's north pole. The points on the sphere corresponding to this icosahedron's 12 vertices are labeled "A" through "L". The north pole is "A" and the south pole is "B". 

![The locations of the 12 original vertices](images/StartingPoints.png "The locations of the 12 original vertices")

The remaining 10 vertices (at latitude +/- arctan(1/2)) are labeled by "peel". A peel is a vertical slice made of 4 faces of the icosahedron. There are 5 peels, each touching the north and south poles.

#TODO diagram of the peels, with vertices labeled on each and world map section in each triangle (might be a bit hard to do but should be possible)
#TODO make clear which points on the edges of a peel belong to it and which don't, and that the north pole and south pole do not belong to any peels

Each vertex has 5 original neighbors, based on which other vertices it is connected to by the edges of the icosahedron.

This is iteration 0 of the lattice. To get the next iteration, we take all edges in the current lattice, bisect them, and draw lines to fill out the new triangular lattice on each local triangle. 

#TODO diagram of this process

New points are created in a particular order. #TODO describe the point ordering

#TODO describe the concept of parent point and how this creates the point's code, and how north and south poles cannot be parent of anything

#TODO describe the concept of directional parent, with diagrams

#TODO describe the concept of watershed and directional watershed, with diagrams

#TODO describe the neighbor direction system, how to navigate among points within a triangular lattice region (no refraction across peel boundary)

#TODO describe how this is complicated by refraction across peel boundary

#TODO describe the arithmetic for getting neighbors

# Working with spherical coordinates

#TODO describe how to get the nearest point code to a latlon within distance tolerance / iteration limit

#TODO describe how to get the latlon of a point code

# Point numbers

#TODO describe birth number and lookup number systems, show some tables of how they correspond to point codes, make it very clear to reader which systems are for what and why they exist and why you'd use them (and if there is no such reason then get rid of that system, e.g. birth order might not be useful)

