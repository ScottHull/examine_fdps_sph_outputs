from math import pi, sqrt, atan
import numpy as np
from src.identify import ParticleMap, ParticleMapFromFiles
from src import centering

start_time = 0
end_time = 1000
number_processes = 100
path = "/scratch/shull4/outfiles"

pm_start = ParticleMapFromFiles(path=path).read(time=start_time)
pm_end = ParticleMapFromFiles(path=path).read(time=end_time)

imp_velocity = 9.3 * 1000  # m/s
M_earth = 5.972 * 10 ** 24  # kg
R_earth = 6371 * 1000
M_moon = 7.35 * 10 ** 22  # kg
G = 6.67 * 10 ** -11
roche = 2.9 * R_earth
start_target = [p for p in pm_start if p.particle_id < 2]
start_impactor = [p for p in pm_start if p.particle_id >= 2]
disk = [p for p in pm_end if p.label == "DISK"]
escape = [p for p in pm_end if p.label == "ESCAPE"]
planet = [p for p in pm_end if p.label == "PLANET"]

total_tar_mass = sum([p.mass for p in start_target])
total_imp_mass = sum([p.mass for p in start_impactor])
total_mass = total_tar_mass + total_imp_mass
imp_mass_ratio = total_imp_mass / total_mass
com_tar_x, com_tar_y, com_tar_z = centering.center_of_mass(
    x_coords=[p.position_vector[0] for p in start_target],
    y_coords=[p.position_vector[1] for p in start_target],
    z_coords=[p.position_vector[2] for p in start_target],
    masses=[p.mass for p in start_target],
)
com_imp_x, com_imp_y, com_imp_z = centering.center_of_mass(
    x_coords=[p.position_vector[0] for p in start_impactor],
    y_coords=[p.position_vector[1] for p in start_impactor],
    z_coords=[p.position_vector[2] for p in start_impactor],
    masses=[p.mass for p in start_impactor],
)
imp_tar_x_offset = com_imp_x - com_tar_x
imp_tar_y_offset = com_imp_y - com_tar_y
impact_angle = atan(imp_tar_y_offset / imp_tar_x_offset) * (180.0 / pi)
print("IMPACTOR:MASS RATIO: {}".format(imp_mass_ratio))
print("IMPACT ANGLE: {}".format(impact_angle))

disk_mass = sum([p.mass for p in disk])
escape_mass = sum([p.mass for p in escape])
planet_mass = sum([p.mass for p in planet])
disk_ang_momentum = sum([p.angular_momentum for p in disk])
satellite_mass = (1.9 * disk_ang_momentum / sqrt(G * M_earth * roche)) - (1.15 * disk_mass) - (1.9 * escape_mass)
print("DISK MASS: {} ({} M_L)\nESCAPE MASS: {} ({} M_L)\nPROTOEARTH MASS: {} ({} M_E)".format(
    disk_mass, disk_mass / M_moon, escape_mass, escape_mass / M_moon, planet_mass, planet_mass / M_earth
))
print("PREDICTED SATELLITE MASS: {} kg ({} lunar masses)".format(satellite_mass, satellite_mass / M_moon))

target_velocity = sum([np.linalg.norm([vx, vy, vz]) for vx, vy, vz in zip(
    [p.velocity_vector[0] for p in start_target],
    [p.velocity_vector[1] for p in start_target],
    [p.velocity_vector[2] for p in start_target]
)])
impactor_velocity = sum([np.linalg.norm([vx, vy, vz]) for vx, vy, vz in zip(
    [p.velocity_vector[0] for p in start_impactor],
    [p.velocity_vector[1] for p in start_impactor],
    [p.velocity_vector[2] for p in start_impactor]
)])
v_imp_verify_tar = target_velocity / ((total_imp_mass / total_mass) * imp_velocity)
v_imp_verify_imp = impactor_velocity / ((total_tar_mass / total_mass) * imp_velocity)
print("V_IMP_TAR: {} (target: {})\nV_IMP_IMP: {} (target: {})".format(v_imp_verify_tar, imp_velocity, v_imp_verify_imp,
                                                                      imp_velocity))
