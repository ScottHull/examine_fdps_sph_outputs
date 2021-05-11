import os
import numpy as np
from src.identify import ParticleMap
from src.combine import CombineFile
from src.report import make_report
import shutil

start_time = 0
end_time = 5000
interval = 50
number_processes = 100

path_to_outputs = "/scratch/shull4/GI"
outfile_plot_path = "/scratch/shull4/outfiles"

paths = [outfile_plot_path]
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
    make_report(particles=particle_map, time=time, to_directory=outfile_plot_path)
