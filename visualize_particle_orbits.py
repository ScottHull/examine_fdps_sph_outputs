from orbitalelements import Particle
from math import pi, sqrt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Visualize:

    def __init__(self, output_file, center=True):
        self.df = pd.read_csv(output_file, sep='\t', skiprows=2, header=None)
        self.particle_id = self.df[1]
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
            self.center = self.__find_center(resolution=1e4, delta_x=1e7, delta_y=1e7, delta_z=1e7)
            self.x = self.x - self.center[0]
            self.y = self.y - self.center[1]
            self.z = self.z - self.center[2]
        self.a = (12713.6 / 2.0) * 1000.0  # present-day equitorial radius of the Earth in m
        self.b = (12756.2 / 2.0) * 1000.0  # present-day polar radius of the Earth in m
        self.__avg_density = 5.5 * 1000

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

    def __earth_ellipsoid_coords(self):
        phi = np.linspace(0, 2 * np.pi, 256).reshape(256, 1)  # the angle of the projection in the xy-plane
        theta = np.linspace(0, np.pi, 256).reshape(-1, 256)  # the angle from the polar axis, ie the polar angle
        x = self.a * np.sin(theta) * np.cos(phi)
        y = self.a * np.sin(theta) * np.sin(phi)
        z = self.b * np.cos(theta)

        return x, y, z

    def __get_points_within_spheroid(self):
        return [p for p in self.coords if
                        -self.a <= p[0] <= self.a and -self.a <= p[1] <= self.a and -self.b <= p[2] <= self.b]

    def __select_random_particle_outside_spheroid(self):
        for index, p in enumerate(self.coords):
            if p[0] <= -self.a or p[0] >= self.a:
                if p[1] <= -self.a or p[1] >= self.a:
                    if p[2] <= -self.b or p[2] >= self.b:
                        v_x = p[3]
                        v_y = p[4]
                        v_z = p[5]
                        mass = p[6]
                        particle = Particle(
                            position_vector=[p[0], p[1], p[2]],
                            velocity_vector=[v_x, v_y, v_z],
                            mass=mass,
                            mass_grav_body=((4.0 / 3.0) * pi * (self.a**2) * self.b) * self.__avg_density
                        )
                        return particle

    def visualize3D(self):
        ax = Axes3D(plt.figure())
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        earth_ellipsoid = self.__earth_ellipsoid_coords()
        ax.plot_surface(earth_ellipsoid[0], earth_ellipsoid[1], earth_ellipsoid[2], color='blue', alpha=0.6)
        particles_in_earth = self.__get_points_within_spheroid()
        rand_selected_particle = self.__select_random_particle_outside_spheroid()
        ax.scatter([p[0] for p in particles_in_earth], [p[1] for p in particles_in_earth],
                   [p[2] for p in particles_in_earth], color="black", marker="+", alpha=0.2)
        ax.scatter(rand_selected_particle.position_vector[0], rand_selected_particle.position_vector[1],
                   rand_selected_particle.position_vector[2], marker="o", s=80, color="red")
        ax.plot([0, rand_selected_particle.position_vector[0]], [1, rand_selected_particle.position_vector[1]],
                [2, rand_selected_particle.position_vector[2]], linewidth=4, color='red', linestyle="--")
        return ax


v = Visualize(output_file="merged_800.dat", center=True)
ax = v.visualize3D()
plt.show()

