import numpy as np


def center_of_mass(x_coords, y_coords, z_coords, masses):
    masses = np.array(masses, dtype=np.float32)
    total_mass = float(np.sum(masses))
    x_center = sum([a * b for a, b in zip(x_coords, masses)]) / total_mass
    y_center = sum([a * b for a, b in zip(y_coords, masses)]) / total_mass
    z_center = sum([a * b for a, b in zip(z_coords, masses)]) / total_mass
    return x_center, y_center, z_center


def __assign_particle_to_profile(profile_list, density_list, particle_position):
    max_index = len(profile_list) - 1
    for index, i in enumerate(profile_list):
        if index < max_index:
            if i <= particle_position < profile_list[index + 1]:
                density_list[index] += 1
                break


def find_center(x, y, z, mass, resolution, delta_x, delta_y, delta_z):
    com = center_of_mass(x_coords=x, y_coords=y, z_coords=z, masses=mass)
    print("Finding center...")
    x_profile = list(np.arange(com[0] - delta_x, com[0] + delta_x, resolution))
    y_profile = list(np.arange(com[1] - delta_y, com[1] + delta_y, resolution))
    z_profile = list(np.arange(com[2] - delta_z, com[2] + delta_z, resolution))

    x_density = [0 for i in x_profile]
    y_density = [0 for i in y_profile]
    z_density = [0 for i in z_profile]

    particles = zip(x, y, z)
    for p in particles:
        particle_x_pos = p[0]
        particle_y_pos = p[1]
        particle_z_pos = p[2]
        __assign_particle_to_profile(profile_list=x_profile, density_list=x_density,
                                     particle_position=particle_x_pos)
        __assign_particle_to_profile(profile_list=y_profile, density_list=y_density,
                                     particle_position=particle_y_pos)
        __assign_particle_to_profile(profile_list=z_profile, density_list=z_density,
                                     particle_position=particle_z_pos)

    x_center = x_profile[x_density.index(max(x_density))]
    y_center = y_profile[y_density.index(max(y_density))]
    z_center = z_profile[z_density.index(max(z_density))]

    print("Found center!")
    return x_center, y_center, z_center
