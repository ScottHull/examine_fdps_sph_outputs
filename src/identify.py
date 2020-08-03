from src.orbitalelements import Particle
from src.centering import find_center, center_of_mass
import pandas as pd
from math import pi, sqrt
import numpy as np


class ParticleMap:

    def __init__(self, output_path, center=True, centering_resolution=1e5, centering_delta=1e7,
                 number_expected_bodies=1):
        self.num_bodies = number_expected_bodies
        self.output = pd.read_csv(output_path, skiprows=2, header=None, delimiter="\t")
        self.com = center_of_mass(x_coords=self.output[3], y_coords=self.output[4],
                                  z_coords=self.output[5], masses=self.output[2])
        if center:
            self.earth_center = find_center(x=self.output[3], y=self.output[4],
                                            z=self.output[5], mass=self.output[2], resolution=centering_resolution,
                                            delta_x=centering_delta,
                                            delta_y=centering_delta, delta_z=centering_delta)
        else:
            self.earth_center = self.com
        self.a = (12713.6 / 2.0) * 1000.0  # present-day equitorial radius of the Earth in m
        self.b = (12756.2 / 2.0) * 1000.0  # present-day polar radius of the Earth in m
        self.oblateness = self.calc_oblateness(a=self.a, b=self.b)
        self.mass_protoearth = self.calc_mass_protoearth(a=self.a, b=self.b)

    def calc_oblateness(self, a, b):
        return (a - b) / b

    def refined_oblateness(self, T_star, T_protoearth, K=0.335):
        numerator = (5.0 / 2.0) * ((T_star / T_protoearth) ** 2)
        denominator = 1.0 + ((5.0 / 2.0) - ((15.0 * K) / 4.0)) ** 2
        return numerator / denominator

    def calc_mass_protoearth(self, a, b, rho=5500):
        return ((4.0 / 3.0) * pi * (a ** 2) * b) * rho

    def collect_all_particles(self):
        print("Collecting particles...")
        particles = []
        grav_mass = self.calc_mass_protoearth(a=self.a, b=self.b)
        for row in self.output.index:
            particle_id = int(self.output[1][row])
            position_vector = [self.output[3][row] - self.earth_center[0], self.output[4][row] - self.earth_center[1],
                               self.output[5][row] - self.earth_center[2]]
            velocity_vector = [self.output[6][row], self.output[7][row], self.output[8][row]]
            try:
                p = Particle(
                    particle_id=particle_id,
                    position_vector=position_vector,
                    velocity_vector=velocity_vector,
                    mass=self.output[2][row],
                    mass_grav_body=grav_mass,
                    entropy=self.output[13][row],
                    temperature=self.output[14][row]
                )
                particles.append(p)
            except:
                pass
        print("Collected particles!")
        return particles

    def __solve(self, particles, planet_label):
        iteration = 0
        K = 0.335
        G = 6.674 * 10 ** -11
        CONVERGENCE = False
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
                    p.assigned_body = planet_label
                else:  # the particle's radial position is not within the planet
                    if p.eccentricity < 1.0:  # elliptic orbit, will remain in the disk
                        if abs(
                                p.periapsis) <= self.a:  # the particle's periapsis is < the equatorial radius and will eventually fall on the planet and become part of the planet
                            NUM_PARTICLES_WITH_PERIAPSES_WITHIN_RADIAL_DISTANCE += 1
                            p.label = "PLANET"
                            p.assigned_body = planet_label
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

    def solve(self):
        print("Beginning solution iteration...")
        target_label = "TARGET"
        impactor_label = "IMPACTOR"
        # particles = self.select_random_particles(max_iteration=50000, max_randint=110000)
        particles = self.collect_all_particles()

        print("Solving target...")
        target = self.__solve(particles=particles, planet_label=target_label)
        print("Finished solving target!")
        if self.num_bodies != 1:
            target_removed_particles = [p for p in particles if p.assigned_body != target_label]
            print("Solving impactor...")
            self.com = center_of_mass(
                x_coords=[p.position_vector[0] for p in target_removed_particles],
                y_coords=[p.position_vector[1] for p in target_removed_particles],
                z_coords=[p.position_vector[2] for p in target_removed_particles],
                masses=[p.mass for p in target_removed_particles]
            )
            impactor = self.__solve(particles=target_removed_particles, planet_label=impactor_label)
            print("Finished solving impactor!")

        return particles
