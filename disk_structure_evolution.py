import pandas as pd
import numpy as np
from src.identify import ParticleMap
from src.combine import CombineFile
from src.report import make_report

start_time = 0
end_time = 1000
interval = 50
number_processes = 100

path_to_outputs = "/scratch/shull4/GI2/"

for time in np.arange(start_time, end_time + interval, interval):
    print("Generating report for time: {}".format(time))
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs)
    pm = ParticleMap(output_path=combined_file, center_on_target_iron=True, plot=True)
    particle_map = pm.solve()
    make_report(particles=particle_map, time=time, to_directory="/particle_outputs")
