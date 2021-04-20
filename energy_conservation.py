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

path_to_outputs = "/scratch/shull4/test"

internal_energy = []
kinetic_energy = []
potential_energy = []
total_energy = []

for time in np.arange(start_time, end_time + interval, interval):
    # combine all outputs for a single time iteration into one file
    print(time)
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t")
    internal_energy_time = 0
    kinetic_energy_time = 0
    potential_energy_time = 0
    total_energy_time = 0
    for row in df.index:
        mass = float(df[2][row])
        v_x = float(df[6][row])
        v_y = float(df[7][row])
        v_z = float(df[8][row])
        v_total = sqrt((v_x ** 2) + (v_y ** 2) + (v_z ** 2))
        potential_energy_i = float(df[12][row])
        internal_energy_i = float(df[10][row])
        kinetic_energy_i = (0.5 * mass * (v_total ** 2)) / mass  # units: kJ/kg
        total_energy_i = potential_energy_i + internal_energy_i + kinetic_energy_i

        internal_energy_time += internal_energy_i
        kinetic_energy_time += kinetic_energy_i
        potential_energy_time += abs(potential_energy_i)
        total_energy_time += total_energy_i

    internal_energy.append(internal_energy_time)
    potential_energy.append(potential_energy_time)
    kinetic_energy.append(kinetic_energy_time)
    total_energy.append(total_energy_time)
    os.remove(f)

fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)
ax1.plot(np.arange(start_time, end_time + interval, interval), kinetic_energy, linewidth=2.0, color='black')
ax2.plot(np.arange(start_time, end_time + interval, interval), potential_energy, linewidth=2.0, color='black')
ax3.plot(np.arange(start_time, end_time + interval, interval), internal_energy, linewidth=2.0, color='black')
ax4.plot(np.arange(start_time, end_time + interval, interval), total_energy, linewidth=2.0, color='black')
ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
ax4.set_yscale('log')
ax4.set_xlabel("Time (Iteration)")
ax1.set_ylabel("Kinetic")
ax2.set_ylabel("Potential (abs. val.)")
ax3.set_ylabel("Internal")
ax4.set_ylabel("Total")
ax1.set_title("Energy Conservation")
ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
plt.savefig("energy.png", format='png')
