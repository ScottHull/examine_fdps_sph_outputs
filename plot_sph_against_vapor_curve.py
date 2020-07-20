import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from src.identify import ParticleMap

vapor_file = "src/phase_data/duniteS_vapour_curve.txt"
vapor_df = pd.read_fwf(vapor_file, header=None, skiprows=1)
temperature = vapor_df[0]
density_solid = vapor_df[1]
density_vapor = vapor_df[2]
pressure = vapor_df[3]
entropy_vaporization = vapor_df[5]

sph_file = "merged_800.dat"
sph_df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")
sph_tag = sph_df[1]
sph_entropy = sph_df[13]
sph_temperature = sph_df[14]
sph_density = sph_df[9]
sph_mass = sph_df[2]
sph_x = sph_df[3]
sph_y = sph_df[4]
sph_z = sph_df[5]

pm = ParticleMap(output_path=sph_file)
particle_map = pm.solve()
planet_particles = [p for p in particle_map if p.label == "PLANET"]

# sph_dunite_temperatures = [i for index, i in enumerate(sph_temperature) if sph_tag[index] % 2 == 0]
# sph_dunite_densities = [i for index, i in enumerate(sph_density) if sph_tag[index] % 2 == 0]
# sph_dunite_entropies = [i for index, i in enumerate(sph_entropy) if sph_tag[index] % 2 == 0]

ax = plt.figure().add_subplot(111)
ax.scatter(
    [p.distance / 1000.0 for p in planet_particles if p.particle_id % 2 == 0],
    [p.temperature for p in planet_particles if p.particle_id % 2 == 0],
    color='blue',
    marker="+",
    label="Dunite"
)
ax.scatter(
    [p.distance / 1000.0 for p in planet_particles if p.particle_id % 2 != 0],
    [p.temperature for p in planet_particles if p.particle_id % 2 != 0],
    color='red',
    marker="+",
    label="Iron"
)
ax.set_xlabel("Radial Distance From Earth Center (km)")
ax.set_ylabel("Temperature")
ax.grid()
ax.legend()

ax2 = plt.figure().add_subplot(111)
ax2.scatter(
    [p.position_vector[0] for p in particle_map if p.label == "PLANET"],
    [p.position_vector[1] for p in particle_map if p.label == "PLANET"],
    color="green",
    marker="+",
    label="Planet"
)
ax2.scatter(
    [p.position_vector[0] for p in particle_map if p.label == "DISK"],
    [p.position_vector[1] for p in particle_map if p.label == "DISK"],
    color="blue",
    marker="+",
    label="Disk"
)
ax2.scatter(
    [p.position_vector[0] for p in particle_map if p.label == "ESCAPE"],
    [p.position_vector[1] for p in particle_map if p.label == "ESCAPE"],
    color="red",
    marker="+",
    label="Escape"
)
e = Ellipse(xy=(0, 0), width=pm.a * 2.0, height=pm.b * 2.0, alpha=0.2, color="blue")
ax2.add_artist(e)
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.grid()
ax2.legend()



plt.show()
