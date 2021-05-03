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
    tag = df['']
    x, y, z = df['x'][row], df['y'][row], df['z'][row]

    vel_vec = [vel_x, vel_y, vel_z]
    rel_vel_vec = [rel_vel_x, rel_vel_y, rel_vel_z]
    pos_vec = [x, y, z]


    vel = np.linalg.norm(vel_vec)
    rel_vel = np.linalg.norm(rel_vel_vec)
    pos = np.linalg.norm(pos_vec)
    data.append((label, vel, rel_vel, pos))


