import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path = "/Users/scotthull/Desktop/0.csv"

df = pd.read_csv(path)

target_silicate = []
target_iron = []
impactor_silicate = []
impactor_iron = []

for row in df.index:
    vel_x, vel_y, vel_z = df['v_x_absolute'][row], df['v_y_absolute'][row], df['v_z_absolute'][row]
    rel_vel_x, rel_vel_y, rel_vel_z = df['v_x_relative'][row], df['v_y_relative'][row], df['v_z_relative'][row]
    tag = df['particle_id'][row]
    label = df['label'][row]
    x, y, z = df['x'][row], df['y'][row], df['z'][row]

    vel_vec = [vel_x, vel_y, vel_z]
    rel_vel_vec = [rel_vel_x, rel_vel_y, rel_vel_z]
    pos_vec = [x, y, z]

    vel = np.linalg.norm(vel_vec)
    rel_vel = np.linalg.norm(rel_vel_vec)
    pos = np.linalg.norm(pos_vec)
    
    if tag == 0:
        target_silicate.append([label, vel, rel_vel, pos, [vel_x, vel_y, vel_z], [rel_vel_x, rel_vel_y, rel_vel_z]])
    elif tag == 1:
        target_iron.append([label, vel, rel_vel, pos, [vel_x, vel_y, vel_z], [rel_vel_x, rel_vel_y, rel_vel_z]])
    elif tag == 2:
        impactor_silicate.append([label, vel, rel_vel, pos, [vel_x, vel_y, vel_z], [rel_vel_x, rel_vel_y, rel_vel_z]])
    elif tag ==3:
        impactor_iron.append([label, vel, rel_vel, pos, [vel_x, vel_y, vel_z], [rel_vel_x, rel_vel_y, rel_vel_z]])

total_particles = target_silicate + target_iron + impactor_silicate + impactor_iron
total_abs_velocity_x = sum(i[4][0] for i in total_particles)
total_abs_velocity_y = sum(i[4][1] for i in total_particles)
total_abs_velocity_z = sum(i[4][2] for i in total_particles)
total_rel_velocity_x = sum(i[5][0] for i in total_particles)
total_rel_velocity_y = sum(i[5][1] for i in total_particles)
total_rel_velocity_z = sum(i[5][2] for i in total_particles)
total_abs_velocity = np.linalg.norm([total_abs_velocity_x, total_abs_velocity_y, total_abs_velocity_z])
total_rel_velocity = np.linalg.norm([total_rel_velocity_x, total_rel_velocity_y, total_rel_velocity_z])

mean_vel_target_silicate = sum([p[1] for p in target_silicate]) / len(target_silicate)
mean_vel_target_iron = sum([p[1] for p in target_iron]) / len(target_iron)
mean_vel_impactor_silicate = sum([p[1] for p in impactor_silicate]) / len(impactor_silicate)
mean_vel_impactor_iron = sum([p[1] for p in impactor_iron]) / len(impactor_iron)

max_vel_target_silicate = max([p[1] for p in target_silicate])
max_vel_target_iron = max([p[1] for p in target_iron])
max_vel_impactor_silicate = max([p[1] for p in impactor_silicate])
max_vel_impactor_iron = max([p[1] for p in impactor_iron])

min_vel_target_silicate = min([p[1] for p in target_silicate])
min_vel_target_iron = min([p[1] for p in target_iron])
min_vel_impactor_silicate = min([p[1] for p in impactor_silicate])
min_vel_impactor_iron = min([p[1] for p in impactor_iron])

print(
    "MEAN V TARGET SILICATE: {}\n"
    "MEAN V TARGET IRON: {}\n"
    "MEAN V IMPACTOR SILICATE: {}\n"
    "MEAN V IMPACTOR IRON: {}\n"
    "MAX V TARGET SILICATE: {}\n"
    "MAX V TARGET IRON: {}\n"
    "MAX V IMPACTOR SILICATE: {}\n"
    "MAX V IMPACTOR IRON: {}\n"
    "MIN V TARGET SILICATE: {}\n"
    "MIN V TARGET IRON: {}\n"
    "MIN V IMPACTOR SILICATE: {}\n"
    "MIN V IMPACTOR IRON: {}\n"
    "TOTAL REL VELOCITY: {}\n"
    "TOTAL ABS VELOCITY: {}".format(
        mean_vel_target_silicate, mean_vel_target_iron,
        mean_vel_impactor_silicate, mean_vel_impactor_iron,
        max_vel_target_silicate, max_vel_target_iron,
        max_vel_impactor_silicate, max_vel_impactor_iron,
        min_vel_target_silicate, min_vel_target_iron,
        min_vel_impactor_silicate, min_vel_impactor_iron,
        total_rel_velocity, total_abs_velocity
    )
)


# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.scatter(
#     [i[3] for i in planet],
#     [i[2] for i in planet],
#     marker="+",
#     color='blue',
#     label="PLANET"
# )
# ax.scatter(
#     [i[3] for i in disk],
#     [i[2] for i in disk],
#     marker="+",
#     color='pink',
#     label="DISK"
# )
# ax.scatter(
#     [i[3] for i in escape],
#     [i[2] for i in escape],
#     marker="+",
#     color='red',
#     label="ESCAPE"
# )
# ax.set_xlabel("Radial Distance")
# ax.set_ylabel("Relative Velocity")
# ax.grid()
# ax.legend()
# ax.set_xlim(0, 1e8)
# # plt.savefig("velocities.png", format='png')
# plt.show()
