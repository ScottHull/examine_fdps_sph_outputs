import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path = "/Users/scotthull/Desktop/2000.csv"

df = pd.read_csv(path)

data = []

for row in df.index:
    vel_x, vel_y, vel_z = df['v_x_absolute'][row], df['v_y_absolute'][row], df['v_z_absolute'][row]
    rel_vel_x, rel_vel_y, rel_vel_z = df['v_x_relative'][row], df['v_y_relative'][row], df['v_z_relative'][row]
    label = df['label'][row]
    x, y, z = df['x'][row], df['y'][row], df['z'][row]

    vel_vec = [vel_x, vel_y, vel_z]
    rel_vel_vec = [rel_vel_x, rel_vel_y, rel_vel_z]
    pos_vec = [x, y, z]

    vel = np.linalg.norm(vel_vec)
    rel_vel = np.linalg.norm(rel_vel_vec)
    pos = np.linalg.norm(pos_vec)
    data.append((label, vel, rel_vel, pos))

escape = [p for p in data if p[0] == "ESCAPE"]
disk = [p for p in data if p[0] == "DISK"]
planet = [p for p in data if p[0] == "PLANET"]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(
    [p[1] for p in planet],
    [p[2] for p in planet],
    marker="+",
    color="blue",
    label="PLANET"
)
ax.scatter(
    [p[1] for p in escape],
    [p[2] for p in escape],
    marker="+",
    color="red",
    label="ESCAPE"
)
ax.scatter(
    [p[1] for p in disk],
    [p[2] for p in disk],
    marker="+",
    color="pink",
    label="DISK"
)
ax.set_xlabel("Absolute Velocity")
ax.set_ylabel("Relative Velocity")
ax.grid()
ax.legend()


fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.scatter(
    [p[3] for p in planet],
    [p[2] for p in planet],
    marker="+",
    color="blue",
    label="PLANET"
)
ax2.scatter(
    [p[3] for p in escape],
    [p[2] for p in escape],
    marker="+",
    color="red",
    label="ESCAPE"
)
ax2.scatter(
    [p[3] for p in disk],
    [p[2] for p in disk],
    marker="+",
    color="pink",
    label="DISK"
)
ax2.set_xlabel("Radial Position")
ax2.set_ylabel("Relative Velocity")
ax2.grid()
ax2.legend()

plt.show()
