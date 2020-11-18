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

path_to_outputs_1 = "/scratch/shull4/GI"
path_to_outputs_2 = "/scratch/shull4/GI2"
entropy_plot_path = os.getcwd() + "/entropy"

paths = [entropy_plot_path]
for i in paths:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

rand_particles = []
rand_selected_particles_indices = []
combined_file = CombineFile(num_processes=number_processes, time=end_time, output_path=path_to_outputs_1).combine()
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
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

d = {1: {}, 2:{}}
for i in rand_selected_particles_indices:
    for key in d.keys():
        d[key].update({i: {
            "distance": [],
            "entropy": [],
            "id": [],
            "times": [],
            "density": [],
            "internal_energy": []
        }})
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
    for i in rand_selected_particles_indices:
        p = [p for p in particle_map if p.particle_name == i][0]
        d[1][i]["distance"].append(p.distance / 1000.0)
        d[1][i]["entropy"].append(p.entropy)
        d[1][i]["id"].append(p.particle_name)
        d[1][i]["times"].append(time)
    # times_1 += [time for t in rand_selected_particles_indices]
    # distances_1 += [particle_map[p].distance / 1000.0 for p in rand_selected_particles_indices]
    # entropies_1 += [particle_map[p].entropy for p in rand_selected_particles_indices]
    # ids += [particle_map[i].particle_name for i in rand_selected_particles_indices]

    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs_2).combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
    particle_map = pm.solve()
    os.remove(f)
    for i in rand_selected_particles_indices:
        p = [p for p in particle_map if p.particle_name == i][0]
        d[2][i]["distance"].append(p.distance / 1000.0)
        d[2][i]["entropy"].append(p.entropy)
        d[2][i]["id"].append(p.particle_name)
        d[2][i]["times"].append(time)
    # times_2 += [time for t in rand_selected_particles_indices]
    # distances_2 += [particle_map[p].distance / 1000.0 for p in rand_selected_particles_indices]
    # entropies_2 += [particle_map[p].entropy for p in rand_selected_particles_indices]


def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


# sc = ax.scatter(
#     times_1,
#     entropies_1,
#     c=ids,
# )
# ax2.scatter(
#     times_2,
#     entropies_2,
#     c=ids,
# )
# cbar = plt.colorbar(sc)
# cbar.set_label("Particle ID")
cmap = get_cmap(len(d[1].keys()))
for i in d[1].keys():
    c = cmap(i)
    ax.plot(
        d[1][i]['times'],
        d[1][i]['entropy'],
        c=c
    )
    ax2.plot(
        d[2][i]['times'],
        d[2][i]['entropy'],
        c=c
    )
ax2.set_xlabel("Time Iteration")
ax.set_ylabel("Entropy")
ax2.set_ylabel("Entropy")
# ax.set_xlim(0, 60000)
ax.set_ylim(0, 10000)
ax2.set_ylim(0, 10000)
ax.grid()
ax2.grid()
fig.savefig("multiple_entropy_as_func_of_time.png", format="png")
