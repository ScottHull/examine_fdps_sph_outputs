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
interval = 500
number_processes = 100
num_rand_particles = 10
entropy_lim = 7000

path_to_outputs = "/scratch/shull4/GI"
entropy_plot_path = os.getcwd() + "/entropy"

paths = [entropy_plot_path]
for i in paths:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

rand_selected_particles_indices = []
combined_file = CombineFile(num_processes=number_processes, time=end_time, output_path=path_to_outputs).combine()
f = os.getcwd() + "/merged_{}.dat".format(end_time)
pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
particle_map = pm.solve()
disk_particles = [p for p in particle_map if p.label == "DISK"]
found_particles = 0
while found_particles < num_rand_particles:
    rand_index = randint(0, len(disk_particles) - 1)
    rand_particle = disk_particles[rand_index]
    if rand_particle.entropy >= entropy_lim and rand_particle.particle_name not in rand_selected_particles_indices:
        rand_selected_particles_indices.append(rand_particle.particle_name)
        found_particles += 1

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)

d = {}
for i in rand_selected_particles_indices:
    d.update({i: {
        "distance": [],
        "entropy": [],
        "id": [],
        "times": []
    }})
times = []
distances = []
entropies = []
ids = []

for time in np.arange(start_time, end_time + interval, interval):
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
    particle_map = pm.solve()
    os.remove(f)
    for i in rand_selected_particles_indices:
        p = [p for p in particle_map if p.particle_name == i][0]
        d[i]["distance"].append(p.distance / 1000.0)
        d[i]["entropy"].append(p.entropy)
        d[i]["id"].append(p.particle_name)
        d[i]["times"].append(time)
    # times += [time for t in rand_selected_particles_indices]
    # distances += [particle_map[p].distance / 1000.0 for p in rand_selected_particles_indices]
    # entropies += [particle_map[p].entropy for p in rand_selected_particles_indices]
    # ids += [particle_map[i].particle_name for i in rand_selected_particles_indices]

# sc = ax.scatter(
#     times,
#     entropies,
#     c=ids,
# )
for i in d.keys():
    ax.plot(
        d[i]['times'],
        d[i]['entropy'],
    )

# cbar = plt.colorbar(sc)
# cbar.set_label("Particle ID")
ax.set_xlabel("Time Iteration")
ax.set_ylabel("Entropy")
# ax.set_xlim(0, 60000)
ax.set_ylim(0, 10000)
ax.grid()
fig.savefig("entropy_as_func_of_time.png", format="png")
