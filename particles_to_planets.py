from math import pi, acos, sqrt
import pandas as pd
from random import randint
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt

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
            return semi_major_axis * (1.0 - eccentricity**2)
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
        return semi_latus_rectum / (1.0 + np.linalg.norm(eccentricity))

    def calc_radius_apoapsis(self, semi_latus_rectum, eccentricity):
        return semi_latus_rectum / (1.0 - np.linalg.norm(eccentricity))

    def calc_rotational_period(self, velocity_vector, orbital_radius):
        return (2.0 * pi * orbital_radius) / np.linalg.norm(velocity_vector)


class Particle(ParticleOrbitalElements):

    def __init__(self, equitorial_radius, particle_mass, particle_id, position_x, position_y, position_z, velocity_x, velocity_y, velocity_z, gravitating_mass):
        super().__init__(position_x, position_y, position_z, velocity_x, velocity_y, velocity_z, gravitating_mass)
        self.particle_mass = particle_mass
        self.particle_id = particle_id
        self.label = self.set_label(equitorial_radius=equitorial_radius)

    def set_label(self, equitorial_radius):
        if abs(np.linalg.norm(self.radius_periapsis)) < equitorial_radius:
            return "PLANET"
        else:
            return "DISK"


class MapParticles:

    def __init__(self, output_path):
        self.output = pd.read_csv(output_path, skiprows=2, header=None, delimiter="\t")
        self.com = self.center_of_mass(x_coords=self.output[3], y_coords=self.output[4],
                                       z_coords=self.output[5], masses=self.output[2])
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

    def plot_particles(self):
        ax = plt.figure().add_subplot(111)
        ax.scatter(self.output[3] - self.com[0], self.output[4] - self.com[1], marker="+", color="black")
        ax.scatter(0, 0, marker="o", s=30, color='red')
        e = Ellipse(xy=(0, 0), width=self.a * 2.0, height=self.b * 2.0, alpha=0.8)
        ax.add_artist(e)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.grid()
        # plt.show()
        return ax

    def plot_particles_from_iteration(self, max_iteration, max_randint):
        ax = self.plot_particles()
        iterated_particles = self.iterate_particles_within_planet(max_iteration, max_randint)
        min_distance, closest_particle = self.__closest_particle_to_equitorial_radius(particles=iterated_particles)
        ax.scatter([i.position_vector[0] for i in iterated_particles],
                   [i.position_vector[1] for i in iterated_particles], marker="+", color="green")
        ax.scatter(closest_particle.position_vector[0], closest_particle.position_vector[1], marker='x',
                   color="blue", s=50)
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

    def iterate_particles_within_planet(self, max_iteration, max_randint):
        particles = []
        for i in range(0, max_iteration):
            r = randint(0, max_randint)
            p = Particle(
                particle_mass=float(self.output[2][r]),
                particle_id=self.output[0][r],
                position_x=self.output[3][r] - self.com[0],
                position_y=self.output[4][r] - self.com[1],
                position_z=self.output[5][r] - self.com[2],
                velocity_x=self.output[6][r],
                velocity_y=self.output[7][r],
                velocity_z=self.output[8][r],
                gravitating_mass=self.calc_mass_protoearth(a=self.a, b=self.b),
                equitorial_radius=self.a
            )
            if p.label == "PLANET":
                particles.append(p)
        return particles

    def __gather_particles(self):
        return [
            Particle(
                particle_mass=self.output[2][row],
                particle_id=self.output[0][row],
                position_x=self.output[3][row] - self.com[0],
                position_y=self.output[4][row] - self.com[1],
                position_z=self.output[5][row] - self.com[2],
                velocity_x=self.output[6][row],
                velocity_y=self.output[7][row],
                velocity_z=self.output[8][row],
                gravitating_mass=self.calc_mass_protoearth(a=self.a, b=self.b),
                equitorial_radius=self.a
            )
            for row in self.output.index
        ]

    def solve(self):
        iteration = 1
        CONVERGENCE = False
        # particles = self.__gather_particles()
        particles = self.iterate_particles_within_planet(max_iteration=5000, max_randint=110000)
        min_distance, closest_particle = self.__closest_particle_to_equitorial_radius(particles=particles)
        print(self.a, min_distance)
        # while CONVERGENCE is False:
        #     L_star = self.__angular_momentum_of_subsection()
        #     I = self.__moment_of_inertia_of_ellpise()
        #     self.__minimum_rotational_period_of_subsection(I=I, L_star=L_star)
            # new_f = self.refined_oblateness()







m = MapParticles(output_path="/Users/scotthull/Desktop/merged_800.dat")
m.plot_particles_from_iteration(max_iteration=5000, max_randint=110000)
# m.solve()