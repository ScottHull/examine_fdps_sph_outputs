import os
import pandas as pd
import numpy as np
from src.identify import ParticleMap
from src.combine import CombineFile
from src.report import make_report
import src.plots as plots
import matplotlib.pyplot as plt
import shutil

start_time = 0
end_time = 1000
interval = 50
number_processes = 100

path_to_outputs = "/scratch/shull4/GI2/"
eccentricity_plot_path = os.getcwd() + "/eccentricity"
disk_structure_path = os.getcwd() + "/structure"
disk_structure_eccentricity_path = os.getcwd() + "/eccentricity_structure"

for time in np.arange(start_time, end_time + interval, interval):
    paths = [eccentricity_plot_path, disk_structure_path, disk_structure_eccentricity_path]
    for i in paths:
        if os.path.exists(i):
            shutil.rmtree(i)
        os.mkdir(i)
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False)
    particle_map = pm.solve()
    # make_report(particles=particle_map, time=time)
    os.remove(f)
    ax = plots.plot_eccentricities(particles=particle_map, a=pm.a, b=pm.b)
    ax.savefig(eccentricity_plot_path + "/{}.png".format(time), format='png')
    plt.close()
    ax = plots.scatter_particles(
                x=[p.position_vector[0] for p in particle_map],
                y=[p.position_vector[1] for p in particle_map],
                tags=[p.particle_id for p in particle_map],
                x_label="x",
                y_label="y",
                a=pm.a,
                b=pm.b
            )
    ax.savefig(disk_structure_path + "/{}.png".format(time), format='png')
    plt.close()
    ax = plots.colorcode_orbits(
                particles=particle_map,
                a=pm.a,
                b=pm.b
            )
    ax.savefig(disk_structure_eccentricity_path + "/{}.png".format(time), format='png')
    plt.close()

plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=eccentricity_plot_path,
              filename="eccentricities.mp4")
plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=disk_structure_path,
              filename="structure.mp4")
plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=disk_structure_eccentricity_path,
              filename="structure_eccentricities.mp4")

