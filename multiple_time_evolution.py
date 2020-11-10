import os
import pandas as pd
import numpy as np
from src.identify import ParticleMap
from src.combine import CombineFile
from src.report import make_report
import src.plots as plots
from src.structure import Structure
import matplotlib.pyplot as plt
import shutil
from random import randint

start_time = 0
end_time = 5000
interval = 50
number_processes = 100
num_rand_particles = 50

path_to_outputs_1 = "/scratch/shull4/GI"
path_to_outputs_2 = "/scratch/shull4/GI2"
entropy_plot_path = os.getcwd() + "/entropy"

paths = [entropy_plot_path]
for i in paths:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

rand_particles = []
combined_file = CombineFile(num_processes=number_processes, time=end_time, output_path=path_to_outputs_1).combine()
f = os.getcwd() + "/merged_{}.dat".format(end_time)
pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
particle_map = pm.solve()
disk_particles = [p for p in particle_map if p.label == "DISK"]
rand_selected_particles_indices = [randint(0, len(disk_particles) - 1) for i in range(0, num_rand_particles)]

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

times_1 = []
distances_1 = []
entropies_1 = []
times_2 = []
distances_2 = []
entropies_2 = []
ids = []

for time in np.arange(start_time, end_time + interval, interval):
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs_1).combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
    particle_map = pm.solve()
    os.remove(f)
    times_1 += [time for t in rand_selected_particles_indices]
    distances_1 += [particle_map[p].distance / 1000.0 for p in rand_selected_particles_indices]
    entropies_1 += [particle_map[p].entropy for p in rand_selected_particles_indices]
    ids += [particle_map[i].particle_name for i in rand_selected_particles_indices]

    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs_2).combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
    particle_map = pm.solve()
    os.remove(f)
    times_2 += [time for t in rand_selected_particles_indices]
    distances_2 += [particle_map[p].distance / 1000.0 for p in rand_selected_particles_indices]
    entropies_2 += [particle_map[p].entropy for p in rand_selected_particles_indices]

sc = ax.scatter(
    times_1,
    entropies_1,
    c=ids,
)
ax2.scatter(
    times_2,
    entropies_2,
    c=ids,
)
cbar = plt.colorbar(sc)
cbar.set_label("Particle ID")
ax.set_xlabel("Time Iteration")
ax.set_ylabel("Entropy")
# ax.set_xlim(0, 60000)
ax.set_ylim(0, 10000)
ax.grid()
fig.savefig("entropy_as_func_of_time.png", format="png")