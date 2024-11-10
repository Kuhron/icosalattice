# import matplotlib.pyplot as plt
# import icosalattice.MapCoordinateMath as mcm

# if __name__ == "__main__":
#     print("testing MapCoordinateMath.py")
#     r_size = 2000
#     c_size = 900
#     r = 1100
#     c = 880
#     # lat00, lon00 = 51.5074, -0.1278  # London
#     # lat01, lon01 = 60.1699, 24.9384  # Helsinki
#     # lat10, lon10 = 40.4168, -3.7038  # Madrid
#     # lat11, lon11 = 41.0082, 28.9784  # Istanbul
#     lat00, lon00 = -33.9249,  18.4241  # Cape Town
#     lat01, lon01 = -31.9505, 115.8605  # Perth
#     lat10, lon10 = -54.8019, -68.3030  # Ushuaia
#     lat11, lon11 = -41.2865, 174.7762  # Wellington
#     p = mcm.get_latlon_of_point_on_map(r, c, r_size, c_size,
#         lat00, lon00, lat01, lon01, lat10, lon10, lat11, lon11, deg=True)
#     # print(p)

#     plt.subplot(111, projection="mollweide")
#     p = plt.plot([-1, 1, 1], [-1, -1, 1], "o-")
#     plt.grid(True)

#     plt.show()
