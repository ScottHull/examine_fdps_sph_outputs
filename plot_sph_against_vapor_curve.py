import pandas as pd
import matplotlib.pyplot as plt

vapor_file = "duniteS_vapour_curve.txt"
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

sph_dunite_temperatures = [i for index, i in enumerate(sph_temperature) if sph_tag[index] % 2 == 0]
sph_dunite_densities = [i for index, i in enumerate(sph_density) if sph_tag[index] % 2 == 0]
sph_dunite_entropies = [i for index, i in enumerate(sph_entropy) if sph_tag[index] % 2 == 0]

ax = plt.figure().add_subplot(111)
ax.scatter(sph_dunite_densities, sph_dunite_entropies, color='blue', marker="+", label="Dunite Particles")
ax.plot(density_vapor, entropy_vaporization, linewidth=1.5, color='red', linestyle="--", label="Vapor Curve (Vapor)")
ax.plot(density_solid, entropy_vaporization, linewidth=1.5, color='purple', linestyle="--", label="Vapor Curve (Vapor)")
ax.set_xlabel("Density")
ax.set_ylabel("Entropy")
ax.grid()
ax.legend()

plt.show()