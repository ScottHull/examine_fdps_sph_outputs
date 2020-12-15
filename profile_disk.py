import os
import numpy as np
import matplotlib.pyplot as plt
from src.identify import ParticleMapFromFiles

# this code works reports.py outputs for pre-profiled simulations
start_time = 0
end_time = 5000
interval = 100
path = "/scratch/shull4/GI_outfiles"

disk_masses = []
escaping_masses = []
planet_masses = []
avg_disk_entropies = []
avg_escaping_entropies = []
disk_entropies = []
escaping_entropies = []
planet_entropies = []
disk_end_positions = []
escaping_end_positions = []
for time in np.arange(start_time, end_time + interval, interval):
    print("At iteration: {}".format(time))
    particle_map = ParticleMapFromFiles(path=path).read(time=time)
    disk_particles = [i for i in particle_map if i.label == "DISK"]
    escaping_particles = [i for i in particle_map if i.label == "ESCAPE"]
    planet_particles = [i for i in particle_map if i.label == "PLANET"]
    planet_masses.append(sum([i.mass for i in planet_particles]))
    disk_masses.append(sum([i.mass for i in disk_particles]))
    escaping_masses.append(sum([i.mass for i in escaping_particles]))
    disk_entropy = [i.entropy for i in disk_particles]
    escape_entropy = [i.entropy for i in escaping_particles]
    planet_entropy = [i.entropy for i in planet_particles]
    if len(disk_particles) > 0:
        avg_disk_entropies.append(sum(disk_entropy) / len(disk_particles))
    else:
        avg_disk_entropies.append(0)
    if len(escaping_particles) > 0:
        avg_escaping_entropies.append(sum(escape_entropy) / len(escaping_particles))
    else:
        avg_escaping_entropies.append(0)
    disk_entropies.append(disk_entropy)
    escaping_entropies.append(escape_entropy)
    planet_entropies.append(planet_entropy)
    if time == end_time:
        disk_end_positions = [i.distance for i in disk_particles]
        escaping_end_positions = [i.distance for i in escaping_particles]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(np.arange(start_time, end_time + interval, interval), planet_masses, linewidth=2.0,
        color="green", label="Planet Mass")
ax.plot(np.arange(start_time, end_time + interval, interval), disk_masses, linewidth=2.0,
        color="blue", label="Disk Mass")
ax.plot(np.arange(start_time, end_time + interval, interval), escaping_masses, linewidth=2.0,
        color="red", label="Escaping Mass")
ax.set_xlabel("Time (iteration)")
ax.set_ylabel("Mass")
ax.set_title("Disk and Escaping Mass")
ax.grid()
ax.legend()
plt.savefig("disk_and_escaping_mass.png", format="png")

fig.clear()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(disk_end_positions, disk_entropies[-1], marker="+",
        color="blue", label="Disk Entropy")
ax.scatter(escaping_end_positions, escaping_entropies[-1], marker="+",
        color="red", label="Escaping Entropy")
ax.set_xlabel("Time (iteration)")
ax.set_ylabel("Entropy")
ax.set_title("Disk and Escaping Entropy at Iteration {}".format(end_time))
ax.grid()
ax.legend()
plt.savefig("disk_and_escaping_entropy.png", format="png")

fig.clear()


