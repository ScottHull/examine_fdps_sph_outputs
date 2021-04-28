import os
import pandas as pd
import numpy as np
from src.identify import ParticleMap, ParticleMapFromFiles
from src.combine import CombineFile
from src.report import make_report
import src.plots_new as plots
from src.structure import Structure
import matplotlib.pyplot as plt
import shutil

start_time = 200
end_time = 2000
interval = 200
number_processes = 100

path_to_outputs = "/scratch/shull4/outfiles"
eccentricity_plot_path = os.getcwd() + "/eccentricity"
disk_structure_path = os.getcwd() + "/structure"
disk_structure_eccentricity_path = os.getcwd() + "/eccentricity_structure"
vmf_path = os.getcwd() + "/vmf"

paths = [eccentricity_plot_path, disk_structure_path, disk_structure_eccentricity_path, vmf_path]
for i in paths:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

for time in np.arange(start_time, end_time + interval, interval):
    pm = ParticleMapFromFiles(path=path_to_outputs)
    particle_map = pm.read(time=time)

    fig = plots.colorcode_orbits(
        particles=particle_map,
        center_plot=True,
        z=True
    )
    fig.savefig(disk_structure_eccentricity_path + "/{}.png".format(time), format='png')

plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=eccentricity_plot_path,
              filename="eccentricities.mp4", fps=5)
plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=disk_structure_path,
              filename="structure.mp4", fps=5)
plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=disk_structure_eccentricity_path,
              filename="structure_eccentricities.mp4", fps=5)
plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=vmf_path,
              filename="vmf.mp4", fps=5)
