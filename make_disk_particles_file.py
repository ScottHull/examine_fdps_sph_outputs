import os
import pandas as pd
import csv
from src.identify import ParticleMap
from src.combine import CombineFile

path_to_outputs = "/scratch/shull4/GI2/"
number_processes = 100
time = 4000

combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
sph_file = os.getcwd() + "/merged_{}.dat".format(time)
sph_df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")
pm = ParticleMap(output_path=sph_file, center_on_target_iron=True, plot=True, relative_velocity=True, center_plot=True)
particle_map = pm.solve()
disk_particles_names = [p.particle_name for p in particle_map if p.label == "DISK"]
print("Disk Particle Count: {}".format(len(disk_particles_names)))

outfile_name = os.getcwd() + "/{}_diskparticles.csv".format(time)
if outfile_name in os.listdir(os.getcwd()):
    os.remove(outfile_name)

outfile = open(outfile_name, 'w')

with open(sph_file, 'r') as infile:
    reader = csv.reader(infile, delimiter="\t")
    time = next(reader)
    header = next(reader)
    print(time)
    print(header)
    header_line = (",".join(str(i) for i in header))
    outfile.write(time[0] + "\n")
    outfile.write(header_line + "\n")
    for row in reader:
        if row[0] in disk_particles_names:
            line = (",".join(str(i) for i in row))
            outfile.write(line + "\n")

os.remove(sph_file)
