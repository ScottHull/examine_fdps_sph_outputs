from math import sqrt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src import centering

path = "/Users/scotthull/Desktop/tar.dat"

df = pd.read_csv(path, sep='\t', skiprows=2, header=None)

tags, mass, x, y, z, vx, vy, vz = df[1], df[2], df[3], df[4], df[5], df[6], df[7], df[8]
radial_distance = [sqrt(i**2 + j**2 + k**2) for i, j, k in zip(x, y, z)]
total_velocity = [sqrt(i**2 + j**2 + k**2) for i, j, k in zip(vx, vy, vz)]

com_x, com_y, com_z = centering.center_of_mass(x_coords=x, y_coords=y, z_coords=z, masses=mass, particle_ids=tags)
min_x, max_x = min(x), max(x)
min_y, max_y = min(y), max(y)
min_z, max_z = min(z), max(z)
mid_x = (max_x + min_x) / 2.0
mid_y = (max_y + min_y) / 2.0
mid_z = (max_z + min_z) / 2.0

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    x,
    y,
    marker="+",
    color='blue'
)
ax.plot(
    [min_x, max_x],
    [mid_y, mid_y],
    linewidth=2.0,
    color='black'
)
ax.plot(
    [mid_x, mid_x],
    [min_y, max_y],
    linewidth=2.0,
    color='black'
)
ax.scatter(
    com_x,
    com_y,
    marker="*",
    s=200,
    color='red'
)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid()

fig2 = plt.figure(figsize=(16, 9))
ax2 = fig2.add_subplot(111)
ax2.scatter(
    [i / 1000.0 for i in radial_distance],
    total_velocity,
    color='black',
    marker="+"
)
ax2.set_xlabel("Radial distance (km)")
ax2.set_ylabel("Total Velocity")
ax2.grid()

plt.show()
