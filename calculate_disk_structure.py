import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from src.identify import ParticleMap
from src.structure import Structure

sph_file = "merged_800.dat"
sph_df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")
pm = ParticleMap(output_path=sph_file)
particle_map = pm.solve()
disk_particles = [p for p in particle_map if p.label == "DISK"]

s = Structure(disk_particles=disk_particles)
vmf = s.calc_vapor_mass_fraction()
print(vmf)
