import os
import pandas as pd
import numpy as np
from src.identify import ParticleMap, ParticleMapFromFiles
from src.combine import CombineFile
from src.report import make_report
import src.plots as plots
from src.structure import Structure
import matplotlib.pyplot as plt
import shutil

time = 50
number_processes = 100
path_to_outputs = "/scratch/shull4/gi"

combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
f = os.getcwd() + "/merged_{}.dat".format(time)
pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
particle_map = pm.solve()

fig = plots.plot_eccentricities(
    particles=particle_map,
    # a=pm.a,
    # b=pm.b
    a=1e6,
    b=1e6,
)