from math import pi, acos, sqrt
import pandas as pd
from random import randint
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import sys

class ParticleOrbitalElements:

    """
    The six orbital elements are as follows:
    - e: eccentricity
    - a: semi-major axis
    - i: inclination
    - Ω: Longitude of the ascending node
    - ω: argument of the periapsis
    - ν: true anomaly
    """

    def __init__(self, position_x, position_y, position_z, velocity_x, velocity_y, velocity_z, gravitating_mass):
        self.position_vector = np.array([position_x, position_y, position_z])  # must be Earth-centered
        self.velocity_vector = np.array([velocity_x, velocity_y, velocity_z])  # must be Earth-centered
        self.gravitating_mass = gravitating_mass
        self.angular_momentum_vector = self.calc_angular_momentum(position_vector=self.position_vector,
                                                           velocity_vector=self.velocity_vector)
        self.magnitude_angular_momentum_xy_plane = np.linalg.norm(self.angular_momentum_vector[0:2])
        self.node_vector = self.calc_node_vector(z_axis_vector=np.array([0, 0, 1]), angular_momentum_vector=self.angular_momentum_vector)
        self.eccentricity_vector = self.calc_orbital_eccentricity_vector(position_vector=self.position_vector,
                                                                         velocity_vector=self.velocity_vector,
                                                                         grav_mass=self.gravitating_mass)
        self.mechanical_energy = self.calc_mechanical_energy(velocity_vector=self.velocity_vector,
                                                             position_vector=self.position_vector,
                                                             grav_mass=self.gravitating_mass)
        self.semi_major_axis = self.calc_semi_major_axis(eccentricity=self.eccentricity_vector,
                                                         specific_mechanical_energy=self.mechanical_energy,
                                                         grav_mass=self.gravitating_mass)
        self.semi_latus_rectum = self.calc_semi_latus_rectum(eccentricity=self.eccentricity_vector,
                                                             angular_momentum_vector=self.angular_momentum_vector,
                                                             grav_mass=self.gravitating_mass,
                                                             semi_major_axis=self.semi_major_axis)
        self.inclination = self.calc_inclination(angular_momentum_vector=self.angular_momentum_vector)
        self.longitude_of_descending_node = self.calc_longitude_of_descending_node(node_vector=self.node_vector)
        self.argument_of_periapsis = self.calc_argument_of_periapsis(node_vector=self.node_vector,
                                                                     eccentricity=self.eccentricity_vector)
        self.true_anomaly = self.calc_true_anomaly(eccentricity=self.eccentricity_vector,
                                                   position_vector=self.position_vector)
        self.checks()
        self.radius_periapsis = self.calc_radius_periapsis(semi_latus_rectum=self.semi_latus_rectum,
                                                           eccentricity=self.eccentricity_vector)
        self.radius_apoapsis = self.calc_radius_apoapsis(semi_latus_rectum=self.semi_latus_rectum,
                                                         eccentricity=self.eccentricity_vector)
        self.rotational_period = self.calc_rotational_period(velocity_vector=self.velocity_vector,
                                                             orbital_radius=self.position_vector)


    def calc_angular_momentum(self, position_vector, velocity_vector):
        return np.cross(position_vector, velocity_vector)

    def calc_node_vector(self, z_axis_vector, angular_momentum_vector):
        return np.cross(z_axis_vector, angular_momentum_vector)

    def graviational_parameter(self, mass):
        G = 6.673 * (10**-11)
        return G * mass

    def calc_orbital_eccentricity_vector(self, position_vector, velocity_vector, grav_mass):
        mu = self.graviational_parameter(mass=grav_mass)
        mag_velocity = np.linalg.norm(velocity_vector)
        mag_position = np.linalg.norm(position_vector)
        numerator = ((mag_velocity**2 - (mu / mag_position)) * position_vector) - \
                    (np.dot(position_vector, velocity_vector) * velocity_vector)
        return numerator / mu

    def calc_mechanical_energy(self, velocity_vector, position_vector, grav_mass):
        mu = self.graviational_parameter(mass=grav_mass)
        mag_velocity = np.linalg.norm(velocity_vector)
        mag_position = np.linalg.norm(position_vector)
        return ((mag_velocity**2) / 2.0) - (mu / mag_position)

    def calc_semi_major_axis(self, eccentricity, specific_mechanical_energy, grav_mass):
        if np.linalg.norm(eccentricity) != 1:
            mu = self.graviational_parameter(mass=grav_mass)
            return - mu / (2 * specific_mechanical_energy)
        else:
            return np.inf()

    def calc_semi_latus_rectum(self, eccentricity, angular_momentum_vector, grav_mass, semi_major_axis):
        if np.linalg.norm(eccentricity) != 1:
            return semi_major_axis * (1.0 - np.linalg.norm(eccentricity)**2)
        else:
            mu = self.graviational_parameter(mass=grav_mass)
            return (np.linalg.norm(angular_momentum_vector)**2) / mu

    def calc_inclination(self, angular_momentum_vector):
        return acos(angular_momentum_vector[2] / np.linalg.norm(angular_momentum_vector)) * (180.0 / pi)

    def calc_longitude_of_descending_node(self, node_vector):
        return acos(node_vector[0] / np.linalg.norm(node_vector)) * (180.0 / pi)

    def calc_argument_of_periapsis(self, node_vector, eccentricity):
        return acos(np.dot(node_vector, eccentricity) / (np.linalg.norm(node_vector) *
                                                         np.linalg.norm(eccentricity))) * (180.0 / pi)

    def calc_true_anomaly(self, eccentricity, position_vector):
        return acos(np.dot(eccentricity, position_vector) / (np.linalg.norm(eccentricity) *
                                                             np.linalg.norm(position_vector))) * (180.0 / pi)

    def checks(self):
        if self.node_vector[1] < 0:
            self.longitude_of_descending_node = 360.0 - self.longitude_of_descending_node
        if np.dot(self.position_vector, self.velocity_vector) < 0:
            self.true_anomaly = 360.0 - self.true_anomaly

    def calc_radius_periapsis(self, semi_latus_rectum, eccentricity):
        # return semi_latus_rectum / (1.0 + np.linalg.norm(eccentricity))
        return self.semi_major_axis * (1.0 - np.linalg.norm(self.eccentricity_vector))

    def calc_radius_apoapsis(self, semi_latus_rectum, eccentricity):
        return semi_latus_rectum / (1.0 - np.linalg.norm(eccentricity))

    def calc_rotational_period(self, velocity_vector, orbital_radius):
        return (2.0 * pi * orbital_radius) / np.linalg.norm(velocity_vector)


class Particle(ParticleOrbitalElements):

    def __init__(self, particle_mass, particle_id, position_x, position_y, position_z, velocity_x, velocity_y, velocity_z, gravitating_mass):
        super().__init__(position_x, position_y, position_z, velocity_x, velocity_y, velocity_z, gravitating_mass)
        self.particle_mass = particle_mass
        self.particle_id = particle_id
        self.particle_radial_distance_from_center = self.__radial_distance_of_particle()
        self.particle_linear_velocity = self.__linear_velocity_of_particle()
        self.label = None

    def __radial_distance_of_particle(self):
        return sqrt(self.position_vector[0]**2 + self.position_vector[1]**2 + self.position_vector[2]**2)

    def __linear_velocity_of_particle(self):
        return sqrt(self.velocity_vector[0]**2 + self.velocity_vector[1]**2 + self.velocity_vector[2]**2)

    def __angular_momentum_of_particle_in_z_direction(self):
        return ((self.position_vector[0] * self.velocity_vector[1]) - (self.position_vector[1] * self.velocity_vector[0]))


class MapParticles:

    def __init__(self, output_path, center=True):
        self.output = pd.read_csv(output_path, skiprows=2, header=None, delimiter="\t")
        self.com = self.center_of_mass(x_coords=self.output[3], y_coords=self.output[4],
                                       z_coords=self.output[5], masses=self.output[2])
        if center:
            self.earth_center = self.__find_center(resolution=5e3, delta_x=1e7, delta_y=1e7, delta_z=1e7)
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
        numerator = (5.0 / 2.0) * ((T_star / T_protoearth)**2)
        denominator = 1.0 + ((5.0 / 2.0) - ((15.0 * K) / 4.0))**2
        return numerator / denominator

    def calc_mass_protoearth(self, a, b, rho=5500):
        return ((4.0 / 3.0) * pi * (a**2) * b) * rho

    def point_in_planet(self, equitorial_radius, polar_radius, x, y, z):
        if abs(z) <= polar_radius:
            r_z_earth = equitorial_radius * sqrt(1.0 - ((abs(z) / polar_radius)**2))
            r_xy = sqrt(x**2 + y**2)
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
        ax.scatter(self.output[3] - self.earth_center[0], self.output[4] - self.earth_center[1], marker="+", color="black")
        if particle_sample is not None:
            ax.scatter([i.position_vector[0] for i in particle_sample],
                       [i.position_vector[1] for i in particle_sample], marker="+", color="green")
        ax.scatter(0, 0, marker="o", s=30, color='red')
        e = Ellipse(xy=(0, 0), width=self.a * 2.0, height=self.b * 2.0, alpha=0.8, color="blue")
        ax.add_artist(e)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.grid()
        # plt.show()
        return ax

    def plot_particles_from_iteration(self, max_iteration, max_randint):
        ax = self.plot_particles()
        iterated_particles = self.select_random_particles(max_iteration, max_randint)

        # min_distance, closest_particle = self.__closest_particle_to_equitorial_radius(particles=iterated_particles)
        ax.scatter([i.position_vector[0] for i in iterated_particles],
                   [i.position_vector[1] for i in iterated_particles], marker="+", color="green")
        # ax.scatter(closest_particle.position_vector[0], closest_particle.position_vector[1], marker='x',
        #            color="blue", s=50)
        plt.show()

    def __moment_of_inertia_of_ellpise(self):
        I_x = (1.0 / 4.0) * (pi * self.a * self.b) * (self.b**2)
        I_y = (1.0 / 4.0) * (pi * self.a * self.b) * (self.a**2)
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
        for i in range(0, max_iteration):
            r = randint(0, max_randint)
            p = Particle(
                particle_mass=float(self.output[2][r]),
                particle_id=self.output[0][r],
                position_x=self.output[3][r] - self.earth_center[0],
                position_y=self.output[4][r] - self.earth_center[1],
                position_z=self.output[5][r] - self.earth_center[2],
                velocity_x=self.output[6][r],
                velocity_y=self.output[7][r],
                velocity_z=self.output[8][r],
                gravitating_mass=self.calc_mass_protoearth(a=self.a, b=self.b),
            )
            particles.append(p)

        print("Collected particles!")
        return particles

    def __gather_particles(self):
        return [
            Particle(
                particle_mass=self.output[2][row],
                particle_id=self.output[0][row],
                position_x=self.output[3][row] - self.earth_center[0],
                position_y=self.output[4][row] - self.earth_center[1],
                position_z=self.output[5][row] - self.earth_center[2],
                velocity_x=self.output[6][row],
                velocity_y=self.output[7][row],
                velocity_z=self.output[8][row],
                gravitating_mass=self.calc_mass_protoearth(a=self.a, b=self.b),
            )
            for row in self.output.index
        ]

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
        # particles = self.__gather_particles()
        particles = self.select_random_particles(max_iteration=5000, max_randint=110000)
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
                if p.particle_radial_distance_from_center <= self.a:  # the particle's radial position is inside of the protoplanetary equatorial radius and is part of the planet
                    NUM_PARTICLES_WITHIN_RADIAL_DISTANCE += 1
                    NEW_MASS_PROTOPLANET += p.particle_mass
                    NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET += p.angular_momentum_vector[2]
                    p.label = "PLANET"
                else:  # the particle's radial position is not within the planet
                    if np.linalg.norm(p.eccentricity_vector) < 1.0:  # elliptic orbit, will remain in the disk
                        if p.radius_periapsis <= self.a:  # the particle's periapsis is < the equatorial radius and will eventually fall on the planet and become part of the planet
                            print(p.radius_periapsis, p.semi_major_axis * (1.0 - np.linalg.norm(p.eccentricity_vector)),
                                  p.semi_latus_rectum / (1.0 + np.linalg.norm(p.eccentricity_vector)))
                            NUM_PARTICLES_WITH_PERIAPSES_WITHIN_RADIAL_DISTANCE += 1
                            p.label = "PLANET"
                            NEW_MASS_PROTOPLANET += p.particle_mass
                            NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET += p.angular_momentum_vector[2]
                        else:  # the particle is part of the disk
                            NUM_PARTICLES_IN_DISK += 1
                            p.label = "DISK"
                            NEW_MASS_DISK += p.particle_mass
                            NEW_Z_ANGULAR_MOMENTUM_DISK += p.angular_momentum_vector[2]
                    else:  # parabolic orbit, will escape the disk
                        NUM_PARTICLES_ESCAPING += 1
                        p.label = "ESCAPE"
                        NEW_MASS_ESCAPED += p.particle_mass
                        NEW_Z_ANGULAR_MOMENTUM_ESCAPED += p.angular_momentum_vector[2]

            moment_of_inertia_protoplanet = (2.0 / 5.0) * NEW_MASS_PROTOPLANET * (self.a**2)
            angular_velocity_protoplanet = NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET / moment_of_inertia_protoplanet
            keplerian_velocity_protoplanet = sqrt((G * NEW_MASS_PROTOPLANET) / self.a**3)
            rotational_period_minimum = (2.0 * pi) / angular_velocity_protoplanet
            rotational_period_protoplanet = (2.0 * pi) / keplerian_velocity_protoplanet
            numerator = (5.0 / 2.0) * ((rotational_period_minimum / rotational_period_protoplanet) ** 2)
            denominator = 1.0 + ((5.0 / 2.0) - ((15.0 * K) / 4.0)) ** 2
            NEW_F = numerator / denominator
            DENSITY_AVG = 5.5 * 1000
            # NEW_A = - self.b / (NEW_F - 1.0)
            NEW_A = ((3 * pi * NEW_MASS_PROTOPLANET * (1 - NEW_F)) / (4 * DENSITY_AVG))**(1/3)
            if (NEW_A - self.a) / self.a < 10**-8:
                CONVERGENCE = True
            else:
                CONVERGENCE = False
            self.a = NEW_A
            self.mass_protoearth = NEW_MASS_PROTOPLANET
            iteration += 1
            print(iteration, (NEW_A - self.a) / self.a)
            print(NEW_MASS_PROTOPLANET, NEW_F, NEW_A, angular_velocity_protoplanet, keplerian_velocity_protoplanet)
            print(
                "NUM_PARTICLES_WITHIN_RADIAL_DISTANCE: {}\n"
                "NUM_PARTICLES_WITH_PERIAPSES_WITHIN_RADIAL_DISTANCE: {}\n"
                "NUM_PARTICLES_IN_DISK: {}\n"
                "NUM_PARTICLES_ESCAPING: {}".format(NUM_PARTICLES_WITHIN_RADIAL_DISTANCE,
                                                    NUM_PARTICLES_WITH_PERIAPSES_WITHIN_RADIAL_DISTANCE,
                                                    NUM_PARTICLES_IN_DISK, NUM_PARTICLES_ESCAPING)
            )
            if iteration == 1:
                ax = self.plot_particles(particle_sample=particles)
                plt.show()
                sys.exit(0)











m = MapParticles(output_path="merged_800.dat", center=True)
# m.plot_particles_from_iteration(max_iteration=5000, max_randint=110000)
m.solve()
# p = m.select_random_particles(max_iteration=5000, max_randint=110000)
# m.find_center(particles=p, resolution=5000, max_x=1e7, max_y=1e7, max_z=1e7, plot_profile=True)