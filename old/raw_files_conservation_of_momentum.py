import os
from math import sqrt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.combine import CombineFile

start_time = 0
end_time = 5000
interval = 200
number_processes = 100
path = "/scratch/shull4/GI_outfiles"


def calc_angular_momentum(velocity_vector, position_vector, mass):
    velocity_vector = [float(velocity_vector[0]), float(velocity_vector[1]), float(velocity_vector[2])]
    position_vector = [float(position_vector[0]), float(position_vector[1]), float(position_vector[2])]
    return float(mass) * np.cross(position_vector, velocity_vector)


angular_momentum_x = []
angular_momentum_y = []
angular_momentum_z = []
angular_momentum_total = []
total_momentum_x = []
total_momentum_y = []
total_momentum_z = []
total_momentums = []
for time in np.arange(start_time, end_time + interval, interval):
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path).combine()
    sph_file = os.getcwd() + "/merged_{}.dat".format(time)
    df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")
    os.remove(sph_file)
    tag, x, y, z, mass, v_x, v_y, v_z = df[1], df[3], df[4], df[5], df[2], df[6], df[7], df[8]
    total_mass = sum([float(m) for m in mass])
    position_vectors = list(zip(x, y, z))
    velocity_vectors = list(zip(v_x, v_y, v_z))
    ams = [calc_angular_momentum(velocity_vectors[index], position_vectors[index], i) for index, i in enumerate(mass)]
    momentum_x = sum([float(m) * float(v) for m, v in zip(mass, v_x)])
    momentum_y = sum([float(m) * float(v) for m, v in zip(mass, v_y)])
    momentum_z = sum([float(m) * float(v) for m, v in zip(mass, v_z)])
    total_momentum = sqrt(momentum_x**2 + momentum_y**2 + momentum_z**2)
    angular_momentum_x.append(sum([i[0] / total_mass for i in ams]))
    angular_momentum_y.append(sum([i[1] / total_mass for i in ams]))
    angular_momentum_z.append(sum([i[2] / total_mass for i in ams]))
    angular_momentum_total.append(sum([sqrt(i[0]**2 + i[1]**2 + i[2]**2) / total_mass for i in ams]))
    total_momentum_x.append(momentum_x / total_mass)
    total_momentum_y.append(momentum_y / total_mass)
    total_momentum_z.append(momentum_z / total_mass)
    total_momentums.append(total_momentum / total_mass)

fig = plt.figure()
fig.set_size_inches(18.5, 10.5)
ax_x = fig.add_subplot(411)
ax_y = fig.add_subplot(412)
ax_z = fig.add_subplot(413)
ax_total = fig.add_subplot(414)
ax_total.set_xlabel("Time Iteration")
ax_x.set_ylabel("Normalized Angular Momentum")
ax_y.set_ylabel("Normalized Angular Momentum")
ax_z.set_ylabel("Normalized Angular Momentum")
ax_total.set_ylabel("Normalized Angular Momentum")
ax_x.set_title("x angular momentum")
ax_y.set_title("y angular momentum")
ax_z.set_title("z angular momentum")
ax_total.set_title("total angular momentum")
ax_x.grid()
ax_y.grid()
ax_z.grid()
ax_total.grid()
ax_x.plot(
    np.arange(start_time, end_time + interval, interval),
    angular_momentum_x,
    linewidth=2.0,
    label="x"
)
ax_y.plot(
    np.arange(start_time, end_time + interval, interval),
    angular_momentum_y,
    linewidth=2.0,
    label="y"
)
ax_z.plot(
    np.arange(start_time, end_time + interval, interval),
    angular_momentum_z,
    linewidth=2.0,
    label="z"
)
ax_total.plot(
    np.arange(start_time, end_time + interval, interval),
    angular_momentum_total,
    linewidth=2.0,
    label="total"
)
plt.savefig("angular_momentum_components.png", format='png')
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
