import os
import pandas as pd
from src.combine import CombineFile
from src.identify import ParticleMap
import matplotlib.pyplot as plt

time = 2000
number_processes = 100
path_to_outputs = "/scratch/shull4/gi"

combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
f = os.getcwd() + "/merged_{}.dat".format(time)
pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True,
                 center_plot=True).collect_all_particles()

impactor = [i for i in pm if i.particle_id > 1]
target = [i for i in pm if i.particle_id <= 1]

min_x_imp = min([i.position_vector[0] for i in impactor])
max_x_imp = max([i.position_vector[0] for i in impactor])
min_y_imp = min([i.position_vector[1] for i in impactor])
max_y_imp = max([i.position_vector[1] for i in impactor])
min_z_imp = min([i.position_vector[2] for i in impactor])
max_z_imp = max([i.position_vector[2] for i in impactor])

mid_x_imp = (max_x_imp + min_x_imp) / 2
mid_y_imp = (max_y_imp + min_y_imp) / 2
mid_z_imp = (max_z_imp + min_z_imp) / 2

min_x_tar = min([i.position_vector[0] for i in target])
max_x_tar = max([i.position_vector[0] for i in target])
min_y_tar = min([i.position_vector[1] for i in target])
max_y_tar = max([i.position_vector[1] for i in target])
min_z_tar = min([i.position_vector[2] for i in target])
max_z_tar = max([i.position_vector[2] for i in target])

mid_x_tar = (max_x_tar + min_x_tar) / 2
mid_y_tar = (max_y_tar + min_y_tar) / 2
mid_z_tar = (max_z_tar + min_z_tar) / 2

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(
    [i.position_vector[0] for i in impactor],
    [i.position_vector[1] for i in impactor],
    marker="+",
    color='red',
    label="IMPACTOR"
)
ax.scatter(
    [i.position_vector[0] for i in target],
    [i.position_vector[1] for i in target],
    marker="+",
    color='red',
    label="TARGET"
)

ax.plot(
    [min_x_imp, max_x_imp],
    [mid_x_imp, mid_x_imp],
    linewidth=2.0,
    color='green'
)
ax.plot(
    [mid_x_imp, mid_x_imp],
    [min_y_imp, max_y_imp],
    linewidth=2.0,
    color='green'
)

ax.plot(
    [min_x_tar, max_x_tar],
    [mid_x_tar, mid_x_tar],
    linewidth=2.0,
    color='pink'
)
ax.plot(
    [mid_x_tar, mid_x_tar],
    [min_y_tar, max_y_tar],
    linewidth=2.0,
    color='pink'
)
ax.set_xlabel("x")
ax.set_xlabel("y")
ax.set_title("Iteration: {}".format(time))
ax.grid()
ax.legend()
plt.savefig("mode1_dimensions.png", format='png')
