from src.interpolation import GenericTrilinearInterpolation
from src.combine import CombineFile
import pandas as pd
import os
import matplotlib.pyplot as plt


path_to_outputs = "/scratch/shull4/impactor/"
number_processes = 100
time = 1000

eos_df_silicate = pd.read_fwf("src/eos/duniteS2.rho_u.txt", header=None, skiprows=2)
density_silicate = list(eos_df_silicate[0])  # load in the full-length density array from eos_df_silicate
energy_silicate = list(eos_df_silicate[1])  # load in the full-length energy array from eos_df_silicate
entropy_silicate = list(eos_df_silicate[5])  # load in the full-length entropy array from eos_df_silicate
tempertaure_silicate = list(eos_df_silicate[2])
pressure_silicate = list(eos_df_silicate[3])
soundspeed_silicate = list(eos_df_silicate[4])

eos_df_iron = pd.read_fwf("src/eos/ironC.rho_u.txt", header=None, skiprows=2)
density_iron = list(eos_df_iron[0])  # load in the full-length density array from eos_df_iron
energy_iron = list(eos_df_iron[1])  # load in the full-length energy array from eos_df_iron
entropy_iron = list(eos_df_iron[5])  # load in the full-length entropy array from eos_df_iron
tempertaure_iron = list(eos_df_iron[2])
pressure_iron = list(eos_df_iron[3])
soundspeed_iron = list(eos_df_iron[4])



combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
sph_file = os.getcwd() + "/merged_{}.dat".format(time)
sph_df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")

given_density = []
given_internal_energy = []
interpolated_internal_energy = []
given_entropy = []
for row in sph_df.index:
    tag_sph = sph_df[1][row]
    density_sph = sph_df[9][row]
    internal_energy_sph = sph_df[10][row]
    entropy_sph = sph_df[13][row]
    given_internal_energy.append(internal_energy_sph)
    given_density.append(density_sph)
    given_entropy.append(entropy_sph)
    if tag_sph % 2 == 0:
        interpolated_val = GenericTrilinearInterpolation(var1_array=density_silicate, var2_array=entropy_silicate,
                                      var3_array=energy_silicate,
                                      var1=density_sph, var2=entropy_sph,
                                      grid_length=120).interpolate()
    else:
        interpolated_val = GenericTrilinearInterpolation(var1_array=density_iron, var2_array=entropy_iron,
                                                         var3_array=energy_iron,
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
