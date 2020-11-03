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
outfile_plot_path = "/scratch/shull4/outfiles"

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
    os.remove(f)

    df = pd.DataFrame({
        "id": [p.particle_name for p in particle_map],
        "tag": [p.particle_id for p in particle_map],
        "distance": [p.distance for p in particle_map],
        "density": [p.density for p in particle_map],
        "internal_energy": [p.internal_energy for p in particle_map],
        "entropy": [p.entropy for p in particle_map],
        "temperature": [p.temperature for p in particle_map],
        "orbit": [p.label for p in particle_map]
    })
    df.to_csv(outfile_plot_path + "/{}.csv".format(time))