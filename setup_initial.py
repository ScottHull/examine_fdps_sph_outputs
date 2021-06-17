from math import pi, sqrt, sin, cos
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src import centering

tar_file = "/Users/scotthull/Desktop/tar.dat"
imp_file = "/Users/scotthull/Desktop/imp.dat"
gi_file = "/Users/scotthull/Desktop/merged_0.dat"
imp_angle = 46.72
imp_angle *= (pi / 180.0)

tar_df = pd.read_csv(tar_file, sep='\t', skiprows=2, header=None)
imp_df = pd.read_csv(imp_file, sep='\t', skiprows=2, header=None)
gi_df = pd.read_csv(gi_file, sep='\t', skiprows=2, header=None)

tar_tags, tar_mass, tar_x, tar_y, tar_z, tar_vx, tar_vy, tar_vz = tar_df[1], tar_df[2], tar_df[3], tar_df[4], tar_df[5], tar_df[6], tar_df[7], tar_df[8]
imp_tags, imp_mass, imp_x, imp_y, imp_z, imp_vx, imp_vy, imp_vz = imp_df[1], imp_df[2], imp_df[3], imp_df[4], imp_df[5], imp_df[6], imp_df[7], imp_df[8]
gi_tags, gi_mass, gi_x, gi_y, gi_z, gi_vx, gi_vy, gi_vz = gi_df[1], gi_df[2], gi_df[3], gi_df[4], gi_df[5], gi_df[6], gi_df[7], gi_df[8]

gi_x_tar = [i for index, i in enumerate(gi_x) if gi_tags[index] < 2]
gi_y_tar = [i for index, i in enumerate(gi_y) if gi_tags[index] < 2]
gi_z_tar = [i for index, i in enumerate(gi_z) if gi_tags[index] < 2]
gi_mass_tar = [i for index, i in enumerate(gi_mass) if gi_tags[index] < 2]
gi_vx_tar = [i for index, i in enumerate(gi_vx) if gi_tags[index] < 2]
gi_vy_tar = [i for index, i in enumerate(gi_vy) if gi_tags[index] < 2]
gi_vz_tar = [i for index, i in enumerate(gi_vz) if gi_tags[index] < 2]
gi_x_imp = [i for index, i in enumerate(gi_x) if gi_tags[index] >= 2]
gi_y_imp = [i for index, i in enumerate(gi_y) if gi_tags[index] >= 2]
gi_z_imp = [i for index, i in enumerate(gi_z) if gi_tags[index] >= 2]
gi_mass_imp = [i for index, i in enumerate(gi_mass) if gi_tags[index] >= 2]
gi_vx_imp = [i for index, i in enumerate(gi_vx) if gi_tags[index] >= 2]
gi_vy_imp = [i for index, i in enumerate(gi_vy) if gi_tags[index] >= 2]
gi_vz_imp = [i for index, i in enumerate(gi_vz) if gi_tags[index] >= 2]

def get_radius_body(coords):
    return max([np.linalg.norm(i) for i in coords])

com_tar_x, com_tar_y, com_tar_z = centering.center_of_mass(
    x_coords=tar_x,
    y_coords=tar_y,
    z_coords=tar_z,
    masses=tar_mass
)
com_imp_x, com_imp_y, com_imp_z = centering.center_of_mass(
    x_coords=imp_x,
    y_coords=imp_y,
    z_coords=imp_z,
    masses=imp_mass
)
# com_tar_gi_x, com_tar_gi_y, com_tar_gi_z = centering.center_of_mass(
#     x_coords=gi_x_tar,
#     y_coords=gi_y_tar,
#     z_coords=gi_z_tar,
#     masses=gi_mass_tar
# )
# com_imp_gi_x, com_imp_gi_y, com_imp_gi_z = centering.center_of_mass(
#     x_coords=gi_x_imp,
#     y_coords=gi_y_imp,
#     z_coords=gi_z_imp,
#     masses=gi_mass_imp
# )

tar_x = [x - com_tar_x for x in tar_x]
tar_y = [y - com_tar_y for y in tar_y]
tar_z = [z - com_tar_z for z in tar_z]
imp_x = [x - com_imp_x for x in imp_x]
imp_y = [y - com_imp_y for y in imp_y]
imp_z = [z - com_imp_z for z in imp_z]
# tar_gi_x = [x - com_tar_gi_x for x in gi_x_tar]
# tar_gi_y = [y - com_tar_gi_y for y in gi_y_tar]
# tar_gi_z = [z - com_tar_gi_z for z in gi_z_tar]
# imp_gi_x = [x - com_imp_gi_x for x in gi_x_imp]
# imp_gi_y = [y - com_imp_gi_y for y in gi_y_imp]
# imp_gi_z = [z - com_imp_gi_z for z in gi_z_imp]

tar_radius = get_radius_body(coords=zip(tar_x, tar_y, tar_z))
imp_radius = get_radius_body(coords=zip(imp_x, imp_y, imp_z))
# gi_tar_radius = get_radius_body(coords=zip(gi_x_tar, gi_y_tar, gi_z_tar))
# gi_imp_radius = get_radius_body(coords=zip(gi_x_imp, gi_y_imp, gi_z_imp))
x_init = cos(imp_angle) * (tar_radius + imp_radius)
y_init = sin(imp_angle) * (tar_radius + imp_radius)
imp_x = [x + x_init for x in imp_x]
imp_y = [y + y_init for y in imp_y]

max_tar_x = max(tar_x)
min_tar_x = min(tar_x)
max_tar_y = max(tar_y)
min_tar_y = min(tar_y)
max_imp_x = max(imp_x)
min_imp_x = min(imp_x)
max_imp_y = max(imp_y)
min_imp_y = min(imp_y)

max_tar_x_gi = max(gi_x_tar)
min_tar_x_gi = min(gi_x_tar)
max_tar_y_gi = max(gi_y_tar)
min_tar_y_gi = min(gi_y_tar)
max_imp_x_gi = max(gi_x_imp)
min_imp_x_gi = min(gi_x_imp)
max_imp_y_gi = max(gi_y_imp)
min_imp_y_gi = min(gi_y_imp)

com_imp_x, com_imp_y, com_imp_z = centering.center_of_mass(
    x_coords=imp_x,
    y_coords=imp_y,
    z_coords=imp_z,
    masses=imp_mass
)

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    tar_x,
    tar_y,
    marker="+",
    color='blue',
    label="target",
    alpha=0.4
)
ax.scatter(
    imp_x,
    imp_y,
    marker="+",
    color='red',
    label="impactor",
    alpha=0.4
)
ax.scatter(
    gi_x_tar,
    gi_y_tar,
    marker="+",
    color='yellow',
    label="target (gi)",
    alpha=0.4
)
ax.scatter(
    gi_x_imp,
    gi_y_imp,
    marker="+",
    color='green',
    label="impactor (gi)",
    alpha=0.4
)
ax.scatter(
    com_tar_x,
    com_tar_y,
    marker="*",
    s=200,
    color='green',
    label="COM TARGET"
)
ax.scatter(
    com_imp_x,
    com_imp_y,
    marker="*",
    s=200,
    color='purple',
    label="COM IMPACTOR"
)
ax.plot(
    [min_tar_x, max_tar_x],
    [(min_tar_y + max_tar_y) / 2, (min_tar_y + max_tar_y) / 2],
    color='red',
    linewidth=2.0,
    linestyle="--"
)
ax.plot(
    [(min_tar_x + max_tar_x) / 2, (min_tar_x + max_tar_x) / 2],
    [min_tar_y, max_tar_y],
    color='red',
    linewidth=2.0,
    linestyle="--"
)
ax.plot(
    [min_imp_x, max_imp_x],
    [(min_imp_y + max_imp_y) / 2, (min_imp_y + max_imp_y) / 2],
    color='red',
    linewidth=2.0,
    linestyle="--"
)
ax.plot(
    [(min_imp_x + max_imp_x) / 2, (min_imp_x + max_imp_x) / 2],
    [min_imp_y, max_imp_y],
    color='red',
    linewidth=2.0,
    linestyle="--"
)
ax.plot(
    [min_tar_x_gi, max_tar_x_gi],
    [(min_tar_y_gi + max_tar_y_gi) / 2, (min_tar_y_gi + max_tar_y_gi) / 2],
    color='blue',
    linewidth=2.0,
    linestyle="--"
)
ax.plot(
    [(min_tar_x_gi + max_tar_x_gi) / 2, (min_tar_x_gi + max_tar_x_gi) / 2],
    [min_tar_y_gi, max_tar_y_gi],
    color='blue',
    linewidth=2.0,
    linestyle="--"
)
ax.plot(
    [min_imp_x_gi, max_imp_x_gi],
    [(min_imp_y_gi + max_imp_y_gi) / 2, (min_imp_y_gi + max_imp_y_gi) / 2],
    color='blue',
    linewidth=2.0,
    linestyle="--"
)
ax.plot(
    [(min_imp_x_gi + max_imp_x_gi) / 2, (min_imp_x_gi + max_imp_x_gi) / 2],
    [min_imp_y_gi, max_imp_y_gi],
    color='blue',
    linewidth=2.0,
    linestyle="--"
)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid()
ax.legend()

plt.show()
