from math import sqrt, acos
import numpy as np

class Particle:

    def __init__(self, position_vector, velocity_vector, mass, mass_grav_body):
        self.__G = 6.674 * 10 ** -11
        self.position_vector = position_vector
        self.velocity_vector = velocity_vector
        self.distance = sqrt(self.position_vector[0]**2 + self.position_vector[1]**2 + self.position_vector[2]**2)
        self.mass = float(mass)
        self.mass_grav_body = mass_grav_body
        self.angular_momentum_vector = self.__angular_momentum()
        self.orbital_energy = self.__total_orbital_energy()
        self.semi_major_axis = self.__semi_major_axis()
        self.eccentricity = self.__eccentricity()
        self.eccentricity_vector = self.__eccentricity_vector()
        self.periapsis_node_vector = self.__node_vector()
        self.inclination = self.__inclination()
        self.longitude_of_ascending_node = self.__longitude_of_ascending_node()
        self.argument_of_periapsis = self.__argument_of_periapsis()
        # self.true_anomaly = self.__true_anomaly()
        self.periapsis = self.__periapsis()


    def __angular_momentum(self):
        return np.cross(self.position_vector, self.velocity_vector)

    def __node_vector(self):
        return np.cross([0, 0, self.position_vector[0]], self.angular_momentum_vector)

    def __total_orbital_energy(self):
        # kinetic energy, KE = 1/2 m v^2
        term1 = (1.0 / 2.0) * self.mass * (np.linalg.norm(self.velocity_vector)**2)
        # vectorized gravitational potential energy, PE = (G M_1 M_2) / r
        term2 = (self.__G * self.mass * self.mass_grav_body) / self.distance
        return term1 - term2

    def __eccentricity(self):
        L = np.linalg.norm(self.angular_momentum_vector)
        E = self.orbital_energy
        alpha = self.__G * self.mass * self.mass_grav_body
        m_reduced = (self.mass * self.mass_grav_body) / (self.mass + self.mass_grav_body)
        return sqrt(1.0 + ((2.0 * E * (L**2)) / (m_reduced * (alpha**2))))

    def __eccentricity_vector(self):
        mu = self.__G * self.mass * self.mass_grav_body
        term1 = ((((np.linalg.norm(self.velocity_vector))**2) / mu) - (1.0 / self.distance)) * \
                np.array(self.position_vector)
        term2 = ((np.dot(self.position_vector, self.velocity_vector) / mu)) * np.array(self.velocity_vector)
        return term1 - term2

    def __semi_major_axis(self):
        E = self.__total_orbital_energy()
        mu = self.__G * self.mass * self.mass_grav_body
        return mu / (2.0 * E)

    def __inclination(self):
        return acos(self.angular_momentum_vector[2] / np.linalg.norm(self.angular_momentum_vector))

    def __longitude_of_ascending_node(self):
        return np.cross([0, 0, 1], self.angular_momentum_vector)

    def __argument_of_periapsis(self):
        return acos(np.dot(self.periapsis_node_vector, self.eccentricity_vector) /
                    (np.linalg.norm(self.periapsis_node_vector) * np.linalg.norm(self.eccentricity_vector)))

    def __true_anomaly(self):
        return acos(np.dot(self.eccentricity_vector, self.position_vector) /
                    (np.linalg.norm(self.eccentricity_vector) * np.linalg.norm(self.position_vector)))

    def __periapsis(self):
        return self.semi_major_axis * (1.0 - self.eccentricity)
