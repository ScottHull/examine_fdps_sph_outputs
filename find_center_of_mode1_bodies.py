import os
from math import pi, sqrt, atan
import pandas as pd
from src.combine import CombineFile
from src.identify import ParticleMap
import matplotlib.pyplot as plt

time = 0
number_processes = 100
path_to_outputs = "/scratch/shull4/gi"

combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
f = os.getcwd() + "/merged_{}.dat".format(time)
pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True,
                 center_plot=True).collect_all_particles(find_orbital_elements=False)


def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
        raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y


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

imp_min_x_coord = (min_x_imp, mid_y_imp)
imp_max_x_coord = (max_x_imp, mid_y_imp)
imp_min_y_coord = (mid_x_imp, min_y_imp)
imp_max_y_coord = (mid_x_imp, max_y_imp)
tar_min_x_coord = (min_x_tar, mid_y_tar)
tar_max_x_coord = (max_x_tar, mid_y_tar)
tar_min_y_coord = (mid_x_tar, min_y_tar)
tar_max_y_coord = (mid_x_tar, max_y_tar)

imp_intersection = line_intersection((imp_min_x_coord, imp_max_x_coord), (imp_min_y_coord, imp_max_y_coord))
tar_intersection = line_intersection((tar_min_x_coord, tar_max_x_coord), (tar_min_y_coord, tar_max_y_coord))

imp_tar_distance = sqrt((imp_intersection[0] - tar_intersection[0])**2 +
                                 (imp_intersection[1] - tar_intersection[1])**2)
imp_tar_x_distance = mid_x_imp - mid_x_tar
imp_tar_y_distance = mid_y_imp - mid_y_tar
imp_tar_angle = atan(imp_tar_y_distance / imp_tar_x_distance) * (180.0 / pi)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(
    [i.position_vector[0] for i in impactor],
    [i.position_vector[1] for i in impactor],
    marker="+",
    color='red',
    label="IMPACTOR",
    alpha=0.8
)
ax.scatter(
    [i.position_vector[0] for i in target],
    [i.position_vector[1] for i in target],
    marker="+",
    color='blue',
    label="TARGET",
    alpha=0.8
)

ax.plot(
    [min_x_imp, max_x_imp],
    [mid_y_imp, mid_y_imp],
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
ax.plot(
    [imp_intersection[0], tar_intersection[0]],
    [imp_intersection[1], tar_intersection[1]],
    linewidth=2.0,
    linestyle="--",
    color='black',
    label="distance: {}\nangle: {} deg".format(imp_tar_distance, imp_tar_angle)
)
ax.plot(
    [mid_x_tar, mid_x_imp],
    [mid_x_tar, mid_x_tar],
    linewidth=2.0,
    linestyle="-.-",
    color='black',
    label="x distance: {}".format(imp_tar_x_distance)
)
ax.plot(
    [mid_x_imp, mid_x_imp],
    [mid_y_tar, mid_y_imp],
    linewidth=2.0,
    linestyle="-.-",
    color='black',
    label="y distance: {}".format(imp_tar_y_distance)
)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Iteration: {}".format(time))
ax.grid()
ax.legend()
plt.savefig("mode1_dimensions.png", format='png')

print(
    "TARGET CENTER COORDS: {}\n"
    "IMPACTOR CENTER COORDS: {}".format(tar_intersection, imp_intersection)
)
