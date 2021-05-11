from orbitalelements import Particle
import numpy as np
import pandas as pd
from random import randint


class Model:

    def __init__(self, output_file, center=True):
        self.df = pd.read_csv(output_file, sep='\t', skiprows=2, header=None)
        self.particle_id = self.df[0]
        self.mass = self.df[2]
        self.x = self.df[3]
        self.y = self.df[4]
        self.z = self.df[5]
        self.v_x = self.df[6]
        self.v_y = self.df[7]
        self.v_z = self.df[8]
        self.coords = list(zip(self.x, self.y, self.z, self.v_x, self.v_y, self.v_z, self.mass))
        self.com = self.center_of_mass(x_coords=self.x, y_coords=self.y, z_coords=self.z, masses=self.mass)
        if center:
            self.center = self.__find_center(resolution=1e5, delta_x=1e7, delta_y=1e7, delta_z=1e7)
            self.x = self.x - self.center[0]
            self.y = self.y - self.center[1]
            self.z = self.z - self.center[2]
        self.a = (12713.6 / 2.0) * 1000.0  # present-day equitorial radius of the Earth in m
        self.b = (12756.2 / 2.0) * 1000.0  # present-day polar radius of the Earth in m
        self.__avg_density = 5.5 * 1000
        self.initial_mass = (4.0 / 3.0) * (self.a ** 2) * self.b * self.__avg_density

    def center_of_mass(self, x_coords, y_coords, z_coords, masses):
        masses = np.array(masses, dtype=np.float32)
        total_mass = float(np.sum(masses))
        x_center = sum([a * b for a, b in zip(x_coords, masses)]) / total_mass
        y_center = sum([a * b for a, b in zip(y_coords, masses)]) / total_mass
        z_center = sum([a * b for a, b in zip(z_coords, masses)]) / total_mass

        return x_center, y_center, z_center

    def __assign_particle_to_profile(self, profile_list, density_list, particle_position):
        max_index = len(profile_list) - 1
        for index, i in enumerate(profile_list):
            if index < max_index:
                if i <= particle_position < profile_list[index + 1]:
                    density_list[index] += 1
                    break

    def __find_center(self, resolution, delta_x, delta_y, delta_z):
        print("Finding center...")

        x_profile = list(np.arange(self.com[0] - delta_x, self.com[0] + delta_x, resolution))
        y_profile = list(np.arange(self.com[1] - delta_y, self.com[1] + delta_y, resolution))
        z_profile = list(np.arange(self.com[2] - delta_z, self.com[2] + delta_z, resolution))

        x_density = [0 for i in x_profile]
        y_density = [0 for i in y_profile]
        z_density = [0 for i in z_profile]

        particles = self.coords
        for p in particles:
            particle_x_pos = p[0]
            particle_y_pos = p[1]
            particle_z_pos = p[2]
            self.__assign_particle_to_profile(profile_list=x_profile, density_list=x_density,
                                              particle_position=particle_x_pos)
            self.__assign_particle_to_profile(profile_list=y_profile, density_list=y_density,
                                              particle_position=particle_y_pos)
            self.__assign_particle_to_profile(profile_list=z_profile, density_list=z_density,
                                              particle_position=particle_z_pos)

        x_center = x_profile[x_density.index(max(x_density))]
        y_center = y_profile[y_density.index(max(y_density))]
        z_center = z_profile[z_density.index(max(z_density))]

        print("Found center!")

        return x_center, y_center, z_center

    def select_random_particle(self):
        len_sample = len(self.x)
        r = randint(0, len_sample)
        rand_coords = [
            self.x[r], self.y[r], self.z[r]
        ]
        rand_velocity = [
            self.v_x[r], self.v_y[r], self.v_z[r]
        ]
        rand_mass = self.mass[r]
        particle_id = self.df[1][r]
        p = Particle(
            position_vector=rand_coords,
            velocity_vector=rand_velocity,
            mass=rand_mass,
            mass_grav_body=self.initial_mass,
            particle_id=particle_id
        )

        return p


m = Model(output_file="merged_800.dat", center=True)
p = m.select_random_particle()
print(p.eccentricity)
print(np.linalg.norm(p.eccentricity_vector))
print(p.semi_major_axis)
