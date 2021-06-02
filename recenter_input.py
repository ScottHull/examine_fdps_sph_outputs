from src import centering
import pandas as pd
import csv
import os

path = "/Users/scotthull/Desktop/tar.dat"

total_N = 0
time = 0
df = pd.read_csv(path, sep='\t', skiprows=2, header=None)

if "input_recenter.csv" in os.listdir(os.getcwd()):
    os.remove("input_recenter.csv")

with open(path, 'r') as infile:
    reader = csv.reader(infile, delimiter="\t")
    time = float(next(reader)[0])
    total_N = int(next(reader)[0])
    infile.close()

tags, mass, x, y, z, vx, vy, vz = df[1], df[2], df[3], df[4], df[5], df[6], df[7], df[8]

com_x, com_y, com_z = centering.center_of_mass(x_coords=x, y_coords=y, z_coords=z, masses=mass, particle_ids=tags)

for row in df.index:
    df[3][row] -= com_x
    df[4][row] -= com_y
    df[5][row] -= com_z
    df[6][row] = 0.0
    df[7][row] = 0.0
    df[8][row] = 0.0

df.to_csv("tmp.dat", index=False, header=False, sep='\t')

outfile = open("input_recenter.csv", 'w')
outfile.write(str(time) + "\n")
outfile.write(str(total_N) + "\n")
with open("tmp.dat", 'r') as infile:
    for line in infile:
        outfile.write(line)
outfile.close()
os.remove("tmp.dat")