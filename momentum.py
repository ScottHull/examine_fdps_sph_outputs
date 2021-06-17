import pandas as pd
import numpy as np

path = "/Users/scotthull/Desktop/merged_0.dat"

df = pd.read_csv(path, sep='\t', skiprows=2, header=None)
tags, mass, x, y, z, vx, vy, vz = df[1], df[2], df[3], df[4], df[5], df[6], df[7], df[8]
norm_velocities = [np.linalg.norm([x, y, z]) for x, y, z in zip(vx, vy, vz)]

target_p = []
impactor_p = []

for index, i in enumerate(tags):
    m = float(mass[index])
    # m = 1
    # p_x, p_y, p_z = m * float(vx[index]), m * float(vy[index]), m * float(vz[index])
    p_x, p_y, p_z = m * float(vx[index]), 0, 0
    total_p = p_x + p_y + p_z
    if i <= 1:
        target_p.append(total_p)
    else:
        impactor_p.append(total_p)

print(
    "TARGET P: {}\nIMPACTOR P: {}\nTOTAL P: {}".format(sum(target_p), sum(impactor_p), sum(target_p) + sum(impactor_p))
)
