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

planet = [i for i in data if i[0] == "PLANET"]
disk = [i for i in data if i[0] == "DISK"]
escape = [i for i in data if i[0] == "ESCAPE"]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(
    [i[3] for i in planet],
    [i[2] for i in planet],
    marker="+",
    color='blue',
    label="PLANET"
)
ax.scatter(
    [i[3] for i in disk],
    [i[2] for i in disk],
    marker="+",
    color='pink',
    label="DISK"
)
ax.scatter(
    [i[3] for i in escape],
    [i[2] for i in escape],
    marker="+",
    color='red',
    label="ESCAPE"
)
ax.set_xlabel("Radial Distance")
ax.set_ylabel("Relative Velocity")
ax.grid()
ax.legend()
ax.set_xlim(0, 1e8)
# plt.savefig("velocities.png", format='png')
plt.show()
