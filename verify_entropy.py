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

start_time = 0
end_time = 5000
interval = 50
number_processes = 100

path_to_outputs = "/scratch/shull4/GI"
entropy_plot_path = os.getcwd() + "/entropy"

paths = [entropy_plot_path]
for i in paths:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

for time in np.arange(start_time, end_time + interval, interval):
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
    particle_map = pm.solve()
    # s = Structure(particles=particle_map, phase="duniteS2")
    # vmf = s.calc_vapor_mass_fraction(target_label="DISK")
    # make_report(particles=particle_map, time=time)
    os.remove(f)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    sc = ax.scatter(
        [p.distance / 1000.0 for p in particle_map],
        [p.entropy for p in particle_map],
        c=[p.particle_id for p in particle_map],
        marker="+"
    )
    cbar = plt.colorbar(sc)
    cbar.set_label("Particle Tag")
    ax.set_xlabel("Distance from Target Center (km)")
    ax.set_ylabel("Entropy")
    ax.set_title("Iteration: {}".format(time))
    ax.set_xlim(0, 20000)
    ax.grid()
    fig.savefig(entropy_plot_path + "/{}.png".format(time), format="png")

plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=entropy_plot_path,
              filename="entropy_evolution.mp4", fps=5)

