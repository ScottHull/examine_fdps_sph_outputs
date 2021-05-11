from src.orbitalelements import Particle
from math import pi, acos, sqrt
import pandas as pd
from random import randint
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import sys


class MapParticles:

    def __init__(self, output_path, center=True, centering_resolution=1e5, centering_delta=1e7):
        self.output = pd.read_csv(output_path, skiprows=2, header=None, delimiter="\t")
        self.com = self.center_of_mass(x_coords=self.output[3], y_coords=self.output[4],
                                       z_coords=self.output[5], masses=self.output[2])
        if center:
            self.earth_center = self.__find_center(resolution=centering_resolution, delta_x=centering_delta,
                                                   delta_y=centering_delta, delta_z=centering_delta)
        else:
            self.earth_center = self.com
        self.a = (12713.6 / 2.0) * 1000.0  # present-day equitorial radius of the Earth in m
        self.b = (12756.2 / 2.0) * 1000.0  # present-day polar radius of the Earth in m
        self.oblateness = self.calc_oblateness(a=self.a, b=self.b)
        self.mass_protoearth = self.calc_mass_protoearth(a=self.a, b=self.b)
        self.tracked_a = [self.a]
        # self.tracked_b = [self.b]
        self.tracked_oblateness = [self.oblateness]
        self.tracked_mass_protoearth = [self.mass_protoearth]

    def calc_oblateness(self, a, b):
        return (a - b) / b

    def refined_oblateness(self, T_star, T_protoearth, K=0.335):
        numerator = (5.0 / 2.0) * ((T_star / T_protoearth) ** 2)
        denominator = 1.0 + ((5.0 / 2.0) - ((15.0 * K) / 4.0)) ** 2
        return numerator / denominator

    def calc_mass_protoearth(self, a, b, rho=5500):
        return ((4.0 / 3.0) * pi * (a ** 2) * b) * rho

    def point_in_planet(self, equitorial_radius, polar_radius, x, y, z):
        if abs(z) <= polar_radius:
            r_z_earth = equitorial_radius * sqrt(1.0 - ((abs(z) / polar_radius) ** 2))
            r_xy = sqrt(x ** 2 + y ** 2)
            if abs(r_xy) <= abs(r_z_earth):
                print(r_xy, r_z_earth)
                return True
        return False

    def center_of_mass(self, x_coords, y_coords, z_coords, masses):
        masses = np.array(masses, dtype=np.float32)
        total_mass = float(np.sum(masses))
        x_center = sum([a * b for a, b in zip(x_coords, masses)]) / total_mass
        y_center = sum([a * b for a, b in zip(y_coords, masses)]) / total_mass
        z_center = sum([a * b for a, b in zip(z_coords, masses)]) / total_mass
        return (x_center, y_center, z_center)

    def plot_particles(self, particle_sample=None):
        ax = plt.figure().add_subplot(111)
        ax.scatter(self.output[3] - self.earth_center[0], self.output[4] - self.earth_center[1], marker="+",
                   color="black")
        if particle_sample is not None:
            ax.scatter([i.position_vector[0] for i in particle_sample if i.label == "PLANET"],
                       [i.position_vector[1] for i in particle_sample if i.label == "PLANET"], marker="+",
                       color="green")
        ax.scatter(0, 0, marker="o", s=30, color='red')
        e = Ellipse(xy=(0, 0), width=self.a * 2.0, height=self.b * 2.0, alpha=0.2, color="blue")
        ax.add_artist(e)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.grid()
        # plt.show()
        return ax

    def __moment_of_inertia_of_ellpise(self):
        I_x = (1.0 / 4.0) * (pi * self.a * self.b) * (self.b ** 2)
        I_y = (1.0 / 4.0) * (pi * self.a * self.b) * (self.a ** 2)
        I = I_x + I_y
        return I

    def __angular_momentum_of_subsection(self):
        K = 0.335
        rho = 5500.0
        G = 6.674 * 10 ** -11
        L_star = K * self.mass_protoearth * self.a * sqrt((G * self.mass_protoearth) / self.a)
        return L_star

    def __minimum_rotational_period_of_subsection(self, I, L_star):
        omega = L_star / I
        T_star = (2.0 * pi) / omega
        return T_star

    def __rotational_period_of_subsection(self, particles):
        return sum([p.magnitude_angular_momentum_xy_plane for p in particles])

    def __closest_particle_to_equitorial_radius(self, particles):
        equitorial_radius = np.linalg.norm([self.a, 0])
        min_distance = None
        min_particle = None
        for index, i in enumerate(particles):
            d = np.linalg.norm([i.position_vector[0], 0])
            if index == 0:
                min_distance = d
                min_particle = i
            else:
                if d < min_distance:
                    min_distance = d
                    min_particle = i
        return min_distance, min_particle

    def __keplerian_velocity(self):
        G = 6.674 * 10 ** -11
        return sqrt((G * self.mass_protoearth / self.a))

    def select_random_particles(self, max_iteration, max_randint):
        print("Collecting particles...")
        particles = []
        grav_mass = self.calc_mass_protoearth(a=self.a, b=self.b)
        for i in range(0, max_iteration):
            r = randint(0, max_randint)
            particle_id = int(self.output[1][r])
            mass = float(self.output[2][r])
            position_vector = [self.output[3][r] - self.earth_center[0], self.output[4][r] - self.earth_center[1],
                               self.output[5][r] - self.earth_center[2]]
            velocity_vector = [self.output[6][r], self.output[7][r], self.output[8][r]]
            p = Particle(
                particle_id=particle_id,
                position_vector=position_vector,
                velocity_vector=velocity_vector,
                mass=mass,
                mass_grav_body=grav_mass
            )
            particles.append(p)
        print("Collected particles!")
        return particles

    def collect_all_particles(self):
        print("Collecting particles...")
        particles = []
        grav_mass = self.calc_mass_protoearth(a=self.a, b=self.b)
        for row in self.output.index:
            particle_id = int(self.output[1][row])
            mass = float(self.output[2][row])
            position_vector = [self.output[3][row] - self.earth_center[0], self.output[4][row] - self.earth_center[1],
                               self.output[5][row] - self.earth_center[2]]
            velocity_vector = [self.output[6][row], self.output[7][row], self.output[8][row]]
            p = Particle(
                particle_id=particle_id,
                position_vector=position_vector,
                velocity_vector=velocity_vector,
                mass=mass,
                mass_grav_body=grav_mass,
                entropy=self.output[13][row],
                temperature=self.output[14][row]
            )
            particles.append(p)
        print("Collected particles!")
        return particles

    def __assign_particle_to_profile(self, profile_list, density_list, particle_position):
        max_index = len(profile_list) - 1
        for index, i in enumerate(profile_list):
            if index < max_index:
                if i <= particle_position < profile_list[index + 1]:
                    density_list[index] += 1
                    break

    def __find_center(self, resolution, delta_x, delta_y, delta_z, plot_profile=False):
        print("Finding center...")
        x_profile = list(np.arange(self.com[0] - delta_x, self.com[0] + delta_x, resolution))
        y_profile = list(np.arange(self.com[1] - delta_y, self.com[1] + delta_y, resolution))
        z_profile = list(np.arange(self.com[2] - delta_z, self.com[2] + delta_z, resolution))

        x_density = [0 for i in x_profile]
        y_density = [0 for i in y_profile]
        z_density = [0 for i in z_profile]

        particles = zip(self.output[3], self.output[4], self.output[5])
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

        if plot_profile:
            ax = plt.figure().add_subplot(111)
            ax.plot(x_profile, x_density, linewidth=2.0, color='black')
            ax.plot(y_profile, y_density, linewidth=2.0, color='red')
            ax.plot(z_profile, z_density, linewidth=2.0, color='blue')
            ax.legend()
            ax.set_xlabel("Position")
            ax.set_ylabel("Number")
            ax.grid()
            plt.show()

        print("Found center!")
        return (x_center, y_center, z_center)

    def solve(self):
        print("Beginning solution iteration...")
        iteration = 0
        CONVERGENCE = False
        K = 0.335
        G = 6.674 * 10 ** -11
        # particles = self.select_random_particles(max_iteration=50000, max_randint=110000)
        particles = self.collect_all_particles()
        while CONVERGENCE is False:
            NUM_PARTICLES_WITHIN_RADIAL_DISTANCE = 0
            NUM_PARTICLES_WITH_PERIAPSES_WITHIN_RADIAL_DISTANCE = 0
            NUM_PARTICLES_IN_DISK = 0
            NUM_PARTICLES_ESCAPING = 0
            NEW_MASS_PROTOPLANET = 0.0
            NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET = 0.0
            NEW_MASS_DISK = 0.0
            NEW_Z_ANGULAR_MOMENTUM_DISK = 0.0
            NEW_MASS_ESCAPED = 0.0
            NEW_Z_ANGULAR_MOMENTUM_ESCAPED = 0.0
            for p in particles:
                if p.distance <= self.a:  # the particle's radial position is inside of the protoplanetary equatorial radius and is part of the planet
                    NUM_PARTICLES_WITHIN_RADIAL_DISTANCE += 1
                    NEW_MASS_PROTOPLANET += p.mass
                    NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET += p.angular_momentum_vector[2]
                    p.label = "PLANET"
                else:  # the particle's radial position is not within the planet
                    if p.eccentricity < 1.0:  # elliptic orbit, will remain in the disk
                        if abs(
                                p.periapsis) <= self.a:  # the particle's periapsis is < the equatorial radius and will eventually fall on the planet and become part of the planet
                            NUM_PARTICLES_WITH_PERIAPSES_WITHIN_RADIAL_DISTANCE += 1
                            p.label = "PLANET"
                            NEW_MASS_PROTOPLANET += p.mass
                            NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET += p.angular_momentum_vector[2]
                        else:  # the particle is part of the disk
                            NUM_PARTICLES_IN_DISK += 1
                            p.label = "DISK"
                            NEW_MASS_DISK += p.mass
                            NEW_Z_ANGULAR_MOMENTUM_DISK += p.angular_momentum_vector[2]
                    else:  # parabolic orbit, will escape the disk
                        NUM_PARTICLES_ESCAPING += 1
                        p.label = "ESCAPE"
                        NEW_MASS_ESCAPED += p.mass
                        NEW_Z_ANGULAR_MOMENTUM_ESCAPED += p.angular_momentum_vector[2]

            moment_of_inertia_protoplanet = (2.0 / 5.0) * NEW_MASS_PROTOPLANET * (self.a ** 2)
            angular_velocity_protoplanet = NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET / moment_of_inertia_protoplanet
            keplerian_velocity_protoplanet = sqrt((G * NEW_MASS_PROTOPLANET) / self.a ** 3)
            rotational_period_minimum = (2.0 * pi) / angular_velocity_protoplanet
            rotational_period_protoplanet = (2.0 * pi) / keplerian_velocity_protoplanet
            numerator = (5.0 / 2.0) * ((angular_velocity_protoplanet / keplerian_velocity_protoplanet) ** 2)
            denominator = 1.0 + ((5.0 / 2.0) - ((15.0 * K) / 4.0)) ** 2
            NEW_F = numerator / denominator
            DENSITY_AVG = 5.5 * 1000
            NEW_A = ((3.0 * pi * NEW_MASS_PROTOPLANET * (1.0 - NEW_F)) / (4.0 * DENSITY_AVG)) ** (1 / 3)
            error = abs((NEW_A - self.a) / self.a)
            if error < 10 ** -8:
                CONVERGENCE = True
            else:
                CONVERGENCE = False
            self.a = NEW_A
            self.mass_protoearth = NEW_MASS_PROTOPLANET
            iteration += 1
            print(
                "ITERATION: {}\n"
                "ERROR: {}\n"
                "NUM_PARTICLES_WITHIN_RADIAL_DISTANCE: {}\n"
                "NUM_PARTICLES_WITH_PERIAPSES_WITHIN_RADIAL_DISTANCE: {}\n"
                "NUM_PARTICLES_IN_DISK: {}\n"
                "NUM_PARTICLES_ESCAPING: {}".format(iteration, error,
                                                    NUM_PARTICLES_WITHIN_RADIAL_DISTANCE,
                                                    NUM_PARTICLES_WITH_PERIAPSES_WITHIN_RADIAL_DISTANCE,
                                                    NUM_PARTICLES_IN_DISK, NUM_PARTICLES_ESCAPING)
            )
            for p in particles:
                try:
                    p.recalculate_elements(mass_grav_body=self.mass_protoearth)
                except:
                    particles.remove(p)
            if CONVERGENCE:
                ax = self.plot_particles(particle_sample=particles)
                plt.show()


m = MapParticles(output_path="merged_1.dat", center=True)
# m.plot_particles_from_iteration(max_iteration=5000, max_randint=110000)
m.solve()
# p = m.select_random_particles(max_iteration=5000, max_randint=110000)
# m.find_center(particles=p, resolution=5000, max_x=1e7, max_y=1e7, max_z=1e7, plot_profile=True)
