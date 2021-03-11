import os
import shutil
from math import sqrt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.combine import CombineFile

"""
This is meant for mode 2 analysis (i.e. single body formation mode).
"""

start_time = 0
end_time = 5000
interval = 50
number_processes = 100

path_to_outputs = "/scratch/shull4/GI"

p_x = []
p_y = []
p_z = []
p_total = []
angular_p = []

for time in np.arange(start_time, end_time + interval, interval):
    # combine all outputs for a single time iteration into one file
    print(time)
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t")
    p_x_time = 0
    p_y_time = 0
    p_z_time = 0
    p_total_time = 0
    angular_p_time = 0
    for row in df.index:
        # iteratve over all particles
        mass = float(df[2][row])
        x = float(df[3][row])
        y = float(df[4][row])
        z = float(df[5][row])
        position_vector = [x, y, z]
        v_x = float(df[6][row])
        v_y = float(df[7][row])
        v_z = float(df[8][row])
        v_total = sqrt((v_x ** 2) + (v_y ** 2) + (v_z ** 2))
        velocity_vector = [v_x, v_y, v_z]
        p_x_i = mass * v_x
        p_y_i = mass * v_y
        p_z_i = mass * v_z
        p_total_i = mass * v_total
        angular_p_i = np.linalg.norm(mass * np.cross(position_vector, velocity_vector))

        # sum particle totals
        p_x_time += p_x_i
        p_y_time += p_y_i
        p_z_time += p_z_i
        p_total_time += p_total_i
        angular_p_time += angular_p_i

    p_x.append(p_x_time)
    p_y.append(p_y_time)
    p_z.append(p_z_time)
    p_total.append(p_total_time)
    angular_p.append(angular_p_time)

    os.remove(f)

fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(np.arange(start_time, end_time + interval, interval), p_total, linewidth=2.0, color='black')
ax1.set_ylabel("Total Momentum (m * v_total)")
ax2.plot(np.arange(start_time, end_time + interval, interval), angular_p, linewidth=2.0, color='black')
ax2.set_ylabel("Angular Momentum (m * [r x v_total])")
ax2.set_xlabel('Time (Iteration)')
ax1.grid()
ax2.grid()
ax1.set_title("Momentum as a Function of Time (Mode 2)")
plt.savefig("momentum.png", format='png')

fig.clear()
fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
ax1.plot(np.arange(start_time, end_time + interval, interval), p_x, linewidth=2.0, color='black')
ax2.plot(np.arange(start_time, end_time + interval, interval), p_y, linewidth=2.0, color='black')
ax3.plot(np.arange(start_time, end_time + interval, interval), p_z, linewidth=2.0, color='black')
ax3.set_xlabel("Time (Iteration)")
ax1.set_ylabel("p_x")
ax2.set_ylabel("p_y")
ax3.set_ylabel("p_z")
ax1.grid()
ax2.grid()
ax3.grid()
ax1.set_title("Momentum Components (Mode 2)")
plt.savefig("momentum_components.png", format='png')
