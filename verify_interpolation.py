from src.combine import CombineFile
from src.interpolation import GenericTrilinearInterpolation
from src.identify import ParticleMap

import os
import pandas as pd
import matplotlib.pyplot as plt

path_to_outputs = "/scratch/shull4/GI2/"
number_processes = 100
time = 5000

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
particles = ParticleMap(output_path=sph_file, center_on_target_iron=False, plot=False, relative_velocity=False,
                        center_plot=False).collect_all_particles()

entropy_interpolation = [
    (p.density, p.internal_energy, p.entropy, GenericTrilinearInterpolation(var1_array=density, var2_array=energy,
                                                                            var3_array=entropy,
                                                                            var1=p.density, var2=p.internal_energy,
                                                                            grid_length=80).interpolate()) for p in
    particles if p.particle_id % 2 == 0]
entropy_errors = [(i[3] - i[2]) / i[2] for i in entropy_interpolation]

for index, i in enumerate(entropy_interpolation):
    if entropy_errors[index] > 0.02
        print([i[0], i[1], i[2], i[3]])

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    [i[0] for index, i in enumerate(entropy_interpolation) if entropy_errors[index] < 100],
    [i[2] for index, i in enumerate(entropy_interpolation) if entropy_errors[index] < 100],
    marker="+",
    color="red",
    label="C++ INTERPOLATION"
)
ax.scatter(
    [i[0] for i in entropy_interpolation],
    [i[3] for i in entropy_interpolation],
    marker="+",
    color="blue",
    label="PYTHON INTERPOLATION"
)
ax.set_xlabel("DENSITY")
ax.set_ylabel("ENTROPY")
ax.set_title("VERIFY ENTROPY INTERPOLATION")
ax.grid()
ax.legend()
fig.savefig("{}_verify_interpolation.png".format(time), format="png")
fig.clear()

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    [i[0] for i in entropy_interpolation],
    [i[3] for i in entropy_interpolation],
    marker="+",
    color="blue",
    label="PYTHON INTERPOLATION"
)
ax.set_xlabel("DENSITY")
ax.set_ylabel("ENTROPY")
ax.set_title("VERIFY PYTHON INTERPOLATION")
ax.grid()
ax.legend()
fig.savefig("{}_python_interpolation.png".format(time), format="png")
fig.clear()

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    [i[0] for index, i in enumerate(entropy_interpolation) if entropy_errors[index] < 100],
    [i[2] for index, i in enumerate(entropy_interpolation) if entropy_errors[index] < 100],
    marker="+",
    color="red",
    label="C++ INTERPOLATION"
)
ax.set_xlabel("DENSITY")
ax.set_ylabel("ENTROPY")
ax.set_title("VERIFY C++ INTERPOLATION")
ax.grid()
ax.legend()
fig.savefig("{}_c_interpolation.png".format(time), format="png")
fig.clear()

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    [i[0] for i in entropy_interpolation],
    [(i[3] - i[2]) / i[2] for i in entropy_interpolation],
    marker="+",
    color="black",
)
ax.set_xlabel("DENSITY")
ax.set_ylabel("PYTHON-C++ INTERPOLATION ERROR")
ax.set_title("PYTHON-C++ INTERPOLATION ERROR")
ax.grid()
fig.savefig("{}_verify_interpolation_error.png".format(time), format="png")
fig.clear()

for i in
