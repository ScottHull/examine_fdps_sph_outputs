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

time = 2000
number_processes = 100
path_to_outputs = "/scratch/shull4/gi"

combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
f = os.getcwd() + "/merged_{}.dat".format(time)
pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
particle_map = pm.solve()
# pm = ParticleMapFromFiles(path=path_to_outputs)
# particle_map = pm.read(time=time)

# fig = plots.plot_eccentricities(
#     particles=particle_map,
#     a=pm.a,
#     b=pm.b
# )
# fig.savefig("eccentricity_{}.png".format(time), format='png')

fig = plots.colorcode_orbits(
        particles=particle_map,
        a=pm.a,
        b=pm.b,
        center_plot=True,
        z=None,
        time=time
    )
fig.savefig("orbits_{}.png".format(time), format='png')
