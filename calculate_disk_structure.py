import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from src.identify import ParticleMap
from src.structure import Structure
from src.combine import CombineFile
import src.plots as plots
import matplotlib.pyplot as plt

path_to_outputs = "/scratch/shull4/GI2/"
number_processes = 100
time = 4000

combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
sph_file = os.getcwd() + "/merged_{}.dat".format(time)
sph_df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")
pm = ParticleMap(output_path=sph_file, center_on_target_iron=True, plot=True, relative_velocity=True, center_plot=True)
particle_map = pm.solve()
disk_particles = [p for p in particle_map if p.label == "DISK"]
os.remove(sph_file)

s = Structure(disk_particles=disk_particles, phase="duniteS2")
vmf = s.calc_vapor_mass_fraction(target_label="DISK")
fig = plots.plot_vfm(
    phase_curve_1_x=s.phase_df['entropy_sol_liq'],
    phase_curve_1_y=s.phase_df['temperature'],
    phase_curve_2_x=s.phase_df['entropy_vap'],
    phase_curve_2_y=s.phase_df['temperature'],
    particles_x=[p.entropy for p in particle_map if p.label == "DISK" and p.particle_id % 2 == 0],
    particles_y=[p.temperature for p in particle_map if p.label == "DISK" and p.particle_id % 2 == 0],
    particles_color=[p.distance for p in particle_map if p.label == "DISK" and p.particle_id % 2 == 0],
    xlabel="Entropy (S)",
    ylabel="Temperature (T [deg K])",
    cbar_label="Radial Distance From Target Center",
    phase_curve_1_label="sol-liq",
    phase_curve_2_label="liq-gas"
)
fig.savefig(os.getcwd() + "/vmf_{}.png".format(time), format="png")
plt.close()
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    [p.distance / 1000.0 for p in particle_map if p.label == "DISK" and p.particle_id % 2 == 0],
    [p.entropy for p in particle_map if p.label == "DISK" and p.particle_id % 2 == 0],
    marker="+",
    color='black'
)
ax.set_xlabel("Radial Distance from Target Center (km)")
ax.set_ylabel("Entropy")
ax.grid()
fig.savefig(os.getcwd() + "/{}_distance_vs_distance.png".format(time), format="png")
print(vmf)
# surface_densities, sorted_distances = s.calc_disk_surface_density()
#
# ax = plt.figure().add_subplot(111)
# ax.plot([i / 1000.0 for i in sorted_distances], surface_densities, color="black", linewidth=2.0)
# ax.set_xlabel("Distance from Earth Center (km)")
# ax.set_ylabel("Disk Surface Density (kg/m3)")
# ax.set_title("Disk + Escaping Particles Surface Density")
# ax.grid()
#
