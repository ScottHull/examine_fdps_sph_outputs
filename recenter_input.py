from src import centering
import pandas as pd

path = "/Users/scotthull/Desktop/imp.dat"
df = pd.read_csv(path, sep='\t', skiprows=2, header=None)

tags, mass, x, y, z, vx, vy, vz = df[1], df[2], df[3], df[4], df[5], df[6], df[7], df[8]

com_x, com_y, com_z = centering.center_of_mass(x_coords=x, y_coords=y, z_coords=z, masses=mass, particle_ids=tags)

for row in df.index:
    df[3][row] -= com_x
    df[4][row] -= com_y
    df[5][row] -= com_z
    df[6][row] = 0.0
    df[7][row] = 0.0
    df[8][row] = 0.0

df.to_csv("input_reform.dat", index=False, header=False, sep='\t')
