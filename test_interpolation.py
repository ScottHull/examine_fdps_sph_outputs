from src.interpolation import GenericTrilinearInterpolation
from src.combine import CombineFile
import pandas as pd
import os
import matplotlib.pyplot as plt


path_to_outputs = "/scratch/shull4/impactor/"
number_processes = 100
time = 1000

eos_df = pd.read_fwf("src/eos/duniteS2.rho_u.txt", header=None, skiprows=2)
density = list(eos_df[0])  # load in the full-length density array from eos_df
energy = list(eos_df[1])  # load in the full-length energy array from eos_df
entropy = list(eos_df[5])  # load in the full-length entropy array from eos_df
tempertaure = list(eos_df[2])
pressure = list(eos_df[3])
soundspeed = list(eos_df[4])

combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
sph_file = os.getcwd() + "/merged_{}.dat".format(time)
sph_df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")

given_density = []
given_internal_energy = []
interpolated_internal_energy = []
given_entropy = []
for row in sph_df.index:
    density_sph = sph_df[9][row]
    internal_energy_sph = sph_df[10][row]
    entropy_sph = sph_df[13][row]
    given_internal_energy.append(internal_energy_sph)
    given_density.append(density_sph)
    given_entropy.append(entropy_sph)
    interpolated_val = GenericTrilinearInterpolation(var1_array=density, var2_array=entropy,
                                  var3_array=energy,
                                  var1=density_sph, var2=entropy_sph,
                                  grid_length=120).interpolate()
    interpolated_internal_energy.append(interpolated_val)


fig = plt.figure()
ax = fig.add_subplot(111)
sc = ax.scatter(
    given_density,
    [(x - y) / y for x, y in zip(interpolated_internal_energy, given_internal_energy)],
    c=given_entropy,
    marker="+"
)
cbar = plt.colorbar(sc)
cbar.set_label("Entropy (from SPH)")
ax.set_xlabel("Density")
ax.set_ylabel("Python-C++ Internal Energy Interpolation Error")
ax.set_title("Internal Energy Interpolation Error")
ax.grid()
plt.savefig("internal_energy_interpolation_error.png", format="png")
