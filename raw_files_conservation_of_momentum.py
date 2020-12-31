import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.combine import CombineFile

start_time = 0
end_time = 5000
interval = 200
number_processes = 100
path = "/scratch/shull4/GI_outfiles"

fig = plt.figure()
fig.set_size_inches(18.5, 10.5)
ax = fig.add_subplot(111)
ax.set_xlabel("Time Iteration")
ax.set_ylabel("Total Angular Momentum of System")
ax.set_title("System Angular Momentum as Function of Time")
ax.grid()

total_ams = []
total_momentum_x = []
total_momentum_y = []
total_momentum_z = []
total_momentums = []
for time in np.arange(start_time, end_time + interval, interval):
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path).combine()
    sph_file = os.getcwd() + "/merged_{}.dat".format(time)
    df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")
    os.remove(sph_file)
    tag, x, y, z, mass, v_x, v_y, v_z = list(df[1]), list(df[3]), list(df[4]), list(df[5]), list(df[2]), list(
        df[6]), list(df[7]), list(df[8])
    total_am = sum(
        [float(m) * np.cross([float(x_coord), float(y_coord), float(z_coord)], [float(vx), float(vy), float(vz)]) for
         m, x_coord, y_coord, z_coord, vx, vy, vz in zip(mass, x, y, z, v_x, v_y, v_z)])
    momentum_x = sum([m * v for m, v in zip(mass, v_x)])
    momentum_y = sum([m * v for m, v in zip(mass, v_y)])
    momentum_z = sum([m * v for m, v in zip(mass, v_z)])
    total_momentum = momentum_x + momentum_y + momentum_z
    total_mass = sum(mass)
    total_ams.append(total_am / total_mass)
    total_momentum_x.append(momentum_x / total_mass)
    total_momentum_y.append(momentum_y / total_mass)
    total_momentum_z.append(momentum_z / total_mass)
    total_momentums.append(total_momentum / total_mass)

ax.plot(
    np.arange(start_time, end_time + interval, interval),
    total_ams,
    linewidth=2.0,
    color="black"
)

plt.savefig("total_angular_momentum.png", format='png')
fig.clear()

fig = plt.figure()
fig.set_size_inches(18.5, 10.5)
ax_x = fig.add_subplot(311)
ax_y = fig.add_subplot(312)
ax_z = fig.add_subplot(313)
ax_z.set_xlabel("Time Iteration")
ax_x.set_ylabel("Normalized Momentum")
ax_y.set_ylabel("Normalized Momentum")
ax_z.set_ylabel("Normalized Momentum")
ax_x.set_title("x-momentum")
ax_y.set_title("y-momentum")
ax_z.set_title("z-momentum")
ax_x.grid()
ax_y.grid()
ax_z.grid()
ax_x.plot(
    np.arange(start_time, end_time + interval, interval),
    total_momentum_x,
    linewidth=2.0,
    label="x"
)
ax_y.plot(
    np.arange(start_time, end_time + interval, interval),
    total_momentum_y,
    linewidth=2.0,
    label="y"
)
ax_z.plot(
    np.arange(start_time, end_time + interval, interval),
    total_momentum_z,
    linewidth=2.0,
    label="z"
)
plt.savefig("total_momentum_components.png", format='png')
fig.clear()

fig = plt.figure()
fig.set_size_inches(18.5, 10.5)
ax = fig.add_subplot(111)
ax.set_xlabel("Time Iteration")
ax.set_ylabel("Total Normalized Momentum of System")
ax.set_title("System Total Momentum as Function of Time")
ax.grid()
ax.plot(
    np.arange(start_time, end_time + interval, interval),
    total_momentums,
    linewidth=2.0,
)
plt.savefig("total_momentum.png", format='png')
fig.clear()
