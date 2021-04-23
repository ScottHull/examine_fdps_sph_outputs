import os
from math import sqrt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.combine import CombineFile
from src.centering import center_of_mass

start_time = 0
end_time = 5000
interval = 1
number_processes = 100
path = "/scratch/shull4/impactor"

coms = []
for time in np.arange(start_time, end_time + interval, interval):
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path).combine()
    sph_file = os.getcwd() + "/merged_{}.dat".format(time)
    df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")
    os.remove(sph_file)
    tag, x, y, z, mass = df[1], df[3], df[4], df[5], df[2]
    com = center_of_mass(x_coords=x, y_coords=y, z_coords=z, masses=mass, particle_ids=tag, target_iron=False)
    coms.append(com)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(
    np.arange(start_time, end_time + interval, interval),
    [sqrt(i[0]**2 + i[1]**2 + i[2]**2) / 1000.0 for i in coms],
    marker="+",
    color="black"
)
ax.set_xlabel("Time Iteration")
ax.set_ylabel("Distance from Original Center (km)")
ax.set_title("Center of Mass Drift in Mode 2")
ax.grid()
plt.savefig("com_drift.png", format="png")
