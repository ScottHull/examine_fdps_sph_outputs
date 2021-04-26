from math import sqrt, acos, sin, cos, exp
import numpy as np


class Particle:

    def __init__(self, position_vector, velocity_vector, mass, density, mass_grav_body, particle_id,
                 temperature, entropy, internal_energy, particle_name=None, relative_velocity_vector=None):
        self.__G = 6.674 * 10 ** -11
        self.label = None
        self.assigned_body = None
        self.particle_name = particle_name
        self.particle_id = particle_id
        self.position_vector = np.array(position_vector)
        self.velocity_vector = np.array(velocity_vector)
        if relative_velocity_vector is None:
            self.relative_velocity_vector = np.array(self.velocity_vector)
        else:
            self.relative_velocity_vector = np.array(relative_velocity_vector)
        self.distance = np.linalg.norm(position_vector)
        self.mass = float(mass)
        self.density = float(density)
        self.temperature = float(temperature)
        self.entropy = float(entropy)
        self.internal_energy = float(internal_energy)
        self.pressure = None
        self.mass_grav_body = mass_grav_body
        self.angular_momentum_vector = self.__angular_momentum()
        self.angular_momentum = 0
        self.momentum_vector = self.__total_momentum_vector()
        self.semi_major_axis = self.__semi_major_axis()
        self.orbital_energy = self.__total_orbital_energy()
        self.eccentricity = self.__eccentricity()
        self.eccentricity_vector = self.__eccentricity_vector()
        self.periapsis_node_vector = self.__node_vector()
        self.inclination = self.__inclination()
        self.longitude_of_ascending_node = self.__longitude_of_ascending_node()
        self.argument_of_periapsis = self.__argument_of_periapsis()
        # self.true_anomaly = self.__true_anomaly()
        self.periapsis = self.__periapsis()
        self.alpha = 0
        self.mass_reduced = 0
        self.kinetic_energy = 0
        self.potential_energy = 0

    def __total_momentum_vector(self):
        m_x = self.mass * self.velocity_vector[0]
        m_y = self.mass * self.velocity_vector[1]
        m_z = self.mass * self.velocity_vector[2]
        return m_x, m_y, m_z

    def __angular_momentum(self):
        am = self.mass * np.cross(self.position_vector, self.relative_velocity_vector)
        self.angular_momentum = sum(am)
        return am

    def __node_vector(self):
        return np.cross([0, 0, self.position_vector[0]], self.angular_momentum_vector)

    def __total_orbital_energy(self):
        # kinetic energy, KE = 1/2 m v^2
        self.kinetic_energy = (1.0 / 2.0) * self.mass * self.relative_velocity_vector ** 2
        # vectorized gravitational potential energy, PE = (G M_1 M_2) / r
        self.potential_energy = - (self.__G * self.mass_grav_body * self.mass) / self.position_vector
        return self.kinetic_energy + self.potential_energy

    def __eccentricity(self):
        L = np.linalg.norm(self.angular_momentum_vector)
        self.alpha = - self.__G * self.mass * self.mass_grav_body
        self.mass_reduced = (self.mass * self.mass_grav_body) / (self.mass + self.mass_grav_body)
        return sqrt(1.0 + ((2.0 * self.orbital_energy * (L ** 2)) / (self.mass_reduced * (self.alpha ** 2))))

    def __eccentricity_vector(self):
        mu = self.__G * self.mass_grav_body
        term1 = ((((np.linalg.norm(self.relative_velocity_vector)) ** 2) / mu) - (1.0 / self.distance)) * \
                np.array(self.position_vector)
        term2 = ((np.dot(self.position_vector, self.relative_velocity_vector) / mu)) * np.array(
            self.relative_velocity_vector)
        return term1 - term2

    def __semi_major_axis(self):
        E_spec = self.__total_orbital_energy() / self.mass
        mu = self.__G * self.mass_grav_body
        return - mu / (2.0 * E_spec)

    def __inclination(self):
        return acos(self.angular_momentum_vector[2] / np.linalg.norm(self.angular_momentum_vector))

    def __longitude_of_ascending_node(self):
        return np.cross([0, 0, 1], self.angular_momentum_vector)

    def __argument_of_periapsis(self):
        return acos(np.dot(self.periapsis_node_vector, self.eccentricity_vector) /
                    (np.linalg.norm(self.periapsis_node_vector) * np.linalg.norm(self.eccentricity_vector)))

    def __true_anomaly(self):
        """
        http://www.braeunig.us/space/orbmech.htm
        The true anomaly is the position of the orbiting body of the ellipse at a given time.
        The true anomaly equals the mean anomaly for a circular orbit.
        :return:
        """
        return acos(np.dot(self.eccentricity_vector, self.position_vector) /
                    (np.linalg.norm(self.eccentricity_vector) * np.linalg.norm(self.position_vector)))

    def __eccentric_anomaly(self):
        return acos(self.position_vector[0] / self.semi_major_axis)

    def __mean_anomaly(self):
        E = self.__eccentric_anomaly()
        return E - (self.eccentricity * sin(E))

    def __periapsis(self):
        return self.semi_major_axis * (1.0 - self.eccentricity)

    def recalculate_elements(self, mass_grav_body):
        self.mass_grav_body = mass_grav_body
        self.angular_momentum_vector = self.__angular_momentum()
        self.momentum_vector = self.__total_momentum_vector()
        self.semi_major_axis = self.__semi_major_axis()
        self.orbital_energy = self.__total_orbital_energy()
        self.eccentricity = self.__eccentricity()
        self.eccentricity_vector = self.__eccentricity_vector()
        self.periapsis_node_vector = self.__node_vector()
        self.inclination = self.__inclination()
        self.longitude_of_ascending_node = self.__longitude_of_ascending_node()
        self.argument_of_periapsis = self.__argument_of_periapsis()
        # self.true_anomaly = self.__true_anomaly()
        self.periapsis = self.__periapsis()
