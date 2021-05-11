import os
import numpy as np
from src.identify import ParticleMap, ParticleMapFromFiles
import src.plots as plots
import matplotlib.pyplot as plt
import shutil

start_time = 0
end_time = 2000
interval = 10
number_processes = 100
path = "/scratch/shull4/outfiles"
plot_path = "/scratch/shull4/map_final_to_all"

if os.path.exists(plot_path):
    shutil.rmtree(plot_path)
os.mkdir(plot_path)

final = {}
pm_end = ParticleMapFromFiles(path=path).read(time=end_time)
for p in pm_end:
    final.update({(p.particle_name, p.particle_id): p})

for time in np.arange(start_time, end_time + interval, interval):
    print("At iteration: {}".format(time))
    pm = ParticleMapFromFiles(path=path).read(time=time)
    planet = [p for p in pm if final[(p.particle_name, p.particle_id)].label == "PLANET"]
    disk = [p for p in pm if final[(p.particle_name, p.particle_id)].label == "DISK"]
    escape = [p for p in pm if final[(p.particle_name, p.particle_id)].label == "ESCAPE"]

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.set_title("Iteration: {}".format(time))
    ax.set_xlim(-1e8, 1e8)
    ax.set_ylim(-1e8, 1e8)
    ax.scatter(
        [p.position_vector[0] for p in planet],
        [p.position_vector[1] for p in planet],
        marker="+",
        color='blue',
        label='PLANET ({})'.format(len(planet))
    )
    ax.scatter(
        [p.position_vector[0] for p in disk],
        [p.position_vector[1] for p in disk],
        marker="+",
        color='pink',
        label='DISK ({})'.format(len(disk))
    )
    ax.scatter(
        [p.position_vector[0] for p in escape],
        [p.position_vector[1] for p in escape],
        marker="+",
        color='red',
        label='ESCAPE ({})'.format(len(escape))
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid()
    ax.legend(loc='upper right')

    fname = plot_path + "/{}.png".format(time)
    plt.savefig(fname, format='png')

plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=plot_path,
              filename="map_final_to_all.mp4")
shutil.rmtree(plot_path)
