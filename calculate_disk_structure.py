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
surface_densities, sorted_distances = s.calc_disk_surface_density()

ax = plt.figure().add_subplot(111)
ax.plot([i / 1000.0 for i in sorted_distances], surface_densities, color="black", linewidth=2.0)
ax.set_xlabel("Distance from Earth Center (km)")
ax.set_ylabel("Disk Surface Density (kg/m3)")
ax.set_title("Disk + Escaping Particles Surface Density")
ax.grid()

plt.show()