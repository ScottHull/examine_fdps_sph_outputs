from math import sqrt
import pandas as pd
import matplotlib.pyplot as plt

target_path = "/Users/scotthull/Desktop/tar 1.dat"
impactor_path = "/Users/scotthull/Desktop/imp 1.dat"

target_df = pd.read_csv(target_path, header=None, skiprows=2, delimiter="\t")
impactor_df = pd.read_csv(impactor_path, header=None, skiprows=2, delimiter="\t")
target_coordinates = zip(target_df[3], target_df[4], target_df[5])
target_distances = [sqrt(i[0]**2 + i[1]**2 + i[2]**2) / 1000.0 for i in target_coordinates]
target_tag = target_df[1]
target_entropy = target_df[13]
target_temperature = target_df[14]
target_internal_energy = target_df[10]
target_pressure = target_df[11]
target_density = target_df[9]
impactor_coordinates = zip(impactor_df[3], impactor_df[4], impactor_df[5])
impactor_distances = [sqrt(i[0]**2 + i[1]**2 + i[2]**2) / 1000.0 for i in impactor_coordinates]
impactor_tag = impactor_df[1]
impactor_entropy = impactor_df[13]
impactor_internal_energy = impactor_df[10]
impactor_pressure = impactor_df[11]
impactor_density = impactor_df[9]
impactor_temperature = impactor_df[14]

target_fig = plt.figure()
target_fig.set_size_inches(18.5, 10.5)
target_ax_density = target_fig.add_subplot(221)
target_ax_pressure = target_fig.add_subplot(222)
target_ax_internal_energy = target_fig.add_subplot(223)
target_ax_entropy = target_fig.add_subplot(224)
target_ax_density.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] % 2 == 0],
    [i for index, i in enumerate(target_density) if target_tag[index] % 2 == 0],
    color="red",
    marker="+",
    label="Silicate"
)
target_ax_density.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] % 2 == 1],
    [i for index, i in enumerate(target_density) if target_tag[index] % 2 == 1],
    color="blue",
    marker="+",
    label="Iron"
)
target_ax_pressure.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] % 2 == 0],
    [i for index, i in enumerate(target_pressure) if target_tag[index] % 2 == 0],
    color="red",
    marker="+",
    label="Silicate"
)
target_ax_pressure.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] % 2 == 1],
    [i for index, i in enumerate(target_pressure) if target_tag[index] % 2 == 1],
    color="blue",
    marker="+",
    label="Iron"
)
target_ax_internal_energy.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] % 2 == 0],
    [i for index, i in enumerate(target_internal_energy) if target_tag[index] % 2 == 0],
    color="red",
    marker="+",
    label="Silicate"
)
target_ax_internal_energy.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] % 2 == 1],
    [i for index, i in enumerate(target_internal_energy) if target_tag[index] % 2 == 1],
    color="blue",
    marker="+",
    label="Iron"
)
target_ax_entropy.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] % 2 == 0],
    [i for index, i in enumerate(target_entropy) if target_tag[index] % 2 == 0],
    color="red",
    marker="+",
    label="Silicate"
)
target_ax_entropy.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] % 2 == 1],
    [i for index, i in enumerate(target_entropy) if target_tag[index] % 2 == 1],
    color="blue",
    marker="+",
    label="Iron"
)
target_ax_density.set_title("Target Density")
target_ax_pressure.set_title("Target Pressure")
target_ax_internal_energy.set_title("Target Internal Energy")
target_ax_entropy.set_title("Target Entropy")
target_ax_internal_energy.set_xlabel("Distance from Center (km)")
target_ax_entropy.set_xlabel("Distance from Center (km)")
target_ax_density.set_ylabel("Density")
target_ax_pressure.set_ylabel("Pressure")
target_ax_internal_energy.set_ylabel("Internal Energy")
target_ax_entropy.set_ylabel("Entropy")
target_ax_density.grid()
target_ax_pressure.grid()
target_ax_internal_energy.grid()
target_ax_entropy.grid()
target_ax_entropy.legend()
plt.tight_layout()

impactor_fig = plt.figure()
impactor_fig.set_size_inches(18.5, 10.5)
impactor_ax_density = impactor_fig.add_subplot(221)
impactor_ax_pressure = impactor_fig.add_subplot(222)
impactor_ax_internal_energy = impactor_fig.add_subplot(223)
impactor_ax_entropy = impactor_fig.add_subplot(224)
impactor_ax_density.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] % 2 == 0],
    [i for index, i in enumerate(impactor_density) if impactor_tag[index] % 2 == 0],
    color="red",
    marker="+",
    label="Silicate"
)
impactor_ax_density.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] % 2 == 1],
    [i for index, i in enumerate(impactor_density) if impactor_tag[index] % 2 == 1],
    color="blue",
    marker="+",
    label="Iron"
)
impactor_ax_pressure.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] % 2 == 0],
    [i for index, i in enumerate(impactor_pressure) if impactor_tag[index] % 2 == 0],
    color="red",
    marker="+",
    label="Silicate"
)
impactor_ax_pressure.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] % 2 == 1],
    [i for index, i in enumerate(impactor_pressure) if impactor_tag[index] % 2 == 1],
    color="blue",
    marker="+",
    label="Iron"
)
impactor_ax_internal_energy.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] % 2 == 0],
    [i for index, i in enumerate(impactor_internal_energy) if impactor_tag[index] % 2 == 0],
    color="red",
    marker="+",
    label="Silicate"
)
impactor_ax_internal_energy.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] % 2 == 1],
    [i for index, i in enumerate(impactor_internal_energy) if impactor_tag[index] % 2 == 1],
    color="blue",
    marker="+",
    label="Iron"
)
impactor_ax_entropy.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] % 2 == 0],
    [i for index, i in enumerate(impactor_entropy) if impactor_tag[index] % 2 == 0],
    color="red",
    marker="+",
    label="Silicate"
)
impactor_ax_entropy.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] % 2 == 1],
    [i for index, i in enumerate(impactor_entropy) if impactor_tag[index] % 2 == 1],
    color="blue",
    marker="+",
    label="Iron"
)
impactor_ax_density.set_title("Impactor Density")
impactor_ax_pressure.set_title("Impactor Pressure")
impactor_ax_internal_energy.set_title("Impactor Internal Energy")
impactor_ax_entropy.set_title("Impactor Entropy")
impactor_ax_internal_energy.set_xlabel("Distance from Center (km)")
impactor_ax_entropy.set_xlabel("Distance from Center (km)")
impactor_ax_density.set_ylabel("Density")
impactor_ax_pressure.set_ylabel("Pressure")
impactor_ax_internal_energy.set_ylabel("Internal Energy")
impactor_ax_entropy.set_ylabel("Entropy")
impactor_ax_density.grid()
impactor_ax_pressure.grid()
impactor_ax_internal_energy.grid()
impactor_ax_entropy.grid()
impactor_ax_entropy.legend()
plt.tight_layout()

temp_fig = plt.figure()
temp_ax_target = temp_fig.add_subplot(121)
temp_ax_impactor = temp_fig.add_subplot(122)
temp_ax_target.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] % 2 == 0],
    [i for index, i in enumerate(target_temperature) if target_tag[index] % 2 == 0],
    color="red",
    marker="+",
    label="Silicate"
)
temp_ax_target.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] % 2 == 1],
    [i for index, i in enumerate(target_temperature) if target_tag[index] % 2 == 1],
    color="blue",
    marker="+",
    label="Iron"
)
temp_ax_impactor.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] % 2 == 0],
    [i for index, i in enumerate(impactor_temperature) if impactor_tag[index] % 2 == 0],
    color="red",
    marker="+",
    label="Silicate"
)
temp_ax_impactor.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] % 2 == 1],
    [i for index, i in enumerate(impactor_temperature) if impactor_tag[index] % 2 == 1],
    color="blue",
    marker="+",
    label="Iron"
)
temp_ax_target.set_title("Target")
temp_ax_impactor.set_title("Impactor")
temp_ax_target.set_xlabel("Distance from Center (km)")
temp_ax_impactor.set_xlabel("Distance from Center (km)")
temp_ax_target.set_ylabel("Temperature")
temp_ax_impactor.set_ylabel("Temperature")
temp_ax_target.grid()
temp_ax_impactor.grid()
temp_ax_impactor.legend()

plt.show()
