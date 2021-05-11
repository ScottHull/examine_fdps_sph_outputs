import matplotlib.pyplot as plt
import pandas as pd
from math import sqrt

target_path = "/Users/scotthull/Desktop/tar.dat"
impactor_path = "/Users/scotthull/Desktop/imp.dat"

target_df = pd.read_csv(target_path, header=None, skiprows=2, delimiter="\t")
impactor_df = pd.read_csv(impactor_path, header=None, skiprows=2, delimiter="\t")
target_coordinates = zip(target_df[3], target_df[4], target_df[5])
target_distances = [sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2) / 1000.0 for i in target_coordinates]
target_tag = target_df[1]
target_entropy = target_df[13]
target_internal_energy = target_df[10]
target_velocities = [sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2) for i in zip(target_df[6], target_df[7], target_df[8])]
impactor_coordinates = zip(impactor_df[3], impactor_df[4], impactor_df[5])
impactor_distances = [sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2) / 1000.0 for i in impactor_coordinates]
impactor_tag = impactor_df[1]
impactor_entropy = impactor_df[13]
impactor_internal_energy = impactor_df[10]
impactor_velocities = [sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2) for i in
                       zip(impactor_df[6], impactor_df[7], impactor_df[8])]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] == 0],
    [i for index, i in enumerate(target_entropy) if target_tag[index] == 0],
    marker="+",
    color="red",
    label="Target Silicate"
)
ax.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] == 1],
    [i for index, i in enumerate(target_entropy) if target_tag[index] == 1],
    marker="+",
    color="blue",
    label="Target Iron"
)
ax.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] == 0],
    [i for index, i in enumerate(impactor_entropy) if impactor_tag[index] == 0],
    marker="+",
    color="green",
    label="Impactor Silicate"
)
ax.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] == 1],
    [i for index, i in enumerate(impactor_entropy) if impactor_tag[index] == 1],
    marker="+",
    color="purple",
    label="Impactor Iron"
)
ax.set_xlabel("Distance from Center (km)")
ax.set_ylabel("entropy")
ax.set_title("Particle Distance vs Internal Energy for Target and Impactor at Equilibrium")
ax.grid()
ax.legend()

fig_v = plt.figure()
ax_target_v = fig_v.add_subplot(121)
ax_impactor_v = fig_v.add_subplot(122)
ax_target_v.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] == 0],
    [i for index, i in enumerate(target_velocities) if target_tag[index] == 0],
    marker="+",
    color="red",
    label="Silicate"
)
ax_target_v.scatter(
    [i for index, i in enumerate(target_distances) if target_tag[index] == 1],
    [i for index, i in enumerate(target_velocities) if target_tag[index] == 1],
    marker="+",
    color="blue",
    label="Iron"
)
ax_impactor_v.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] == 0],
    [i for index, i in enumerate(impactor_velocities) if impactor_tag[index] == 0],
    marker="+",
    color="red",
    label="Silicate"
)
ax_impactor_v.scatter(
    [i for index, i in enumerate(impactor_distances) if impactor_tag[index] == 1],
    [i for index, i in enumerate(impactor_velocities) if impactor_tag[index] == 1],
    marker="+",
    color="blue",
    label="Iron"
)
ax_impactor_v.set_xlabel("Distance from Center (km)")
ax_target_v.set_xlabel("Distance from Center (km)")
ax_impactor_v.set_ylabel("Velocity")
ax_target_v.set_ylabel("Velocity")
ax_target_v.set_title("Target")
ax_impactor_v.set_title("Impactor")
ax_target_v.grid()
ax_impactor_v.grid()
ax_impactor_v.legend()

plt.show()
