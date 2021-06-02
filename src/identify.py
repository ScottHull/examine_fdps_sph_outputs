from src.orbitalelements import Particle
from src.centering import find_center, center_of_mass
import src.plots as plots
import pandas as pd
from math import pi, sqrt
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import statistics

EARTH_MASS = 5.972 * 10 ** 24
LUNAR_MASS = 7.34767309 * 10 ** 22
L_EM = 3.5 * 10 ** 34


class ParticleMap:

    def __init__(self, output_path, center=True, centering_resolution=1e5, centering_delta=1e7,
                 number_expected_bodies=1, center_on_target_iron=False, plot=False, relative_velocity=True,
                 center_plot=False):
        self.centering_resolution = centering_resolution
        self.centering_delta = centering_delta
        self.num_bodies = number_expected_bodies
        self.__relative_velocity = relative_velocity
        self.__center_on_target_iron = center_on_target_iron
        self.__center_plot = center_plot
        self.output = pd.read_csv(output_path, skiprows=2, header=None, delimiter="\t")
        if center_plot or center_on_target_iron:
            self.com = center_of_mass(x_coords=self.output[3], y_coords=self.output[4],
                                      z_coords=self.output[5], masses=self.output[2], particle_ids=self.output[1],
                                      target_iron=self.__center_on_target_iron)
        else:
            self.com = [0, 0, 0]
        if center:
            self.earth_center = find_center(x=self.output[3], y=self.output[4],
                                            z=self.output[5], mass=self.output[2], particle_ids=self.output[1],
                                            resolution=self.centering_resolution,
                                            delta_x=self.centering_delta, delta_y=self.centering_delta,
                                            delta_z=self.centering_delta,
                                            target_iron_centering=self.__center_on_target_iron)
        else:
            self.earth_center = self.com
        self.a = (12713.6 / 2.0) * 1000.0  # present-day equitorial radius of the Earth in m
        self.b = (12756.2 / 2.0) * 1000.0  # present-day polar radius of the Earth in m
        self.oblateness = self.calc_oblateness(a=self.a, b=self.b)
        self.mass_protoearth = self.calc_mass_protoearth(a=self.a, b=self.b)
        self.plot = plot

    def calc_oblateness(self, a, b):
        return (a - b) / b

    def refined_oblateness(self, T_star, T_protoearth, K=0.335):
        numerator = (5.0 / 2.0) * ((T_star / T_protoearth) ** 2)
        denominator = 1.0 + ((5.0 / 2.0) - ((15.0 * K) / 4.0)) ** 2
        return numerator / denominator

    def calc_mass_protoearth(self, a, b, rho=5500):
        return ((4.0 / 3.0) * pi * (a ** 2) * b) * rho

    def collect_all_particles(self, find_orbital_elements=True):
        print("Collecting particles...")
        particles = []
        grav_mass = self.calc_mass_protoearth(a=self.a, b=self.b)
        for row in self.output.index:
            position_vector = [self.output[3][row] - self.earth_center[0], self.output[4][row] - self.earth_center[1],
                               self.output[5][row] - self.earth_center[2]]
            velocity_vector = [self.output[6][row], self.output[7][row], self.output[8][row]]
            p = Particle(
                particle_name=int(self.output[0][row]),
                particle_id=int(self.output[1][row]),
                position_vector=position_vector,
                velocity_vector=velocity_vector,
                mass=self.output[2][row],
                density=self.output[9][row],
                internal_energy=self.output[10][row],
                mass_grav_body=grav_mass,
                entropy=self.output[13][row],
                temperature=self.output[14][row],
                find_orbital_elements=find_orbital_elements,
            )
            p.pressure = self.output[11][row]

            particles.append(p)

        if self.__relative_velocity:
            self.target_velocity = [
                statistics.mean([p.velocity_vector[0] for p in particles if p.particle_id == 1]),
                statistics.mean([p.velocity_vector[1] for p in particles if p.particle_id == 1]),
                statistics.mean([p.velocity_vector[2] for p in particles if p.particle_id == 1])
            ]
            for p in particles:
                velocity_vector = p.velocity_vector

                relative_velocity_vector = [
                    velocity_vector[0] - self.target_velocity[0],
                    velocity_vector[1] - self.target_velocity[1],
                    velocity_vector[2] - self.target_velocity[2]
                ]
                p.relative_velocity_vector = relative_velocity_vector
                p.recalculate_elements(self.mass_protoearth)
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
                if abs(p.distance) <= self.a:
                    # if abs(p.position_vector[0]) <= self.a and abs(p.position_vector[2]) <= self.a and abs(
                    #         p.position_vector[
                    #             1]) <= self.b:  # the particle's radial position is inside of the protoplanetary polar and equatorial radii and is part of the planet
                    NUM_PARTICLES_WITHIN_RADIAL_DISTANCE += 1
                    NEW_MASS_PROTOPLANET += p.mass
                    NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET += p.angular_momentum_vector[2]
                    p.label = "PLANET"
                    p.assigned_body = planet_label
                else:  # the particle's radial position is not within the planet
                    if p.eccentricity <= 1.0:  # elliptic orbit, will remain in the disk
                        if abs(
                                p.periapsis) < self.a:  # the particle's periapsis is < the equatorial radius and will eventually fall on the planet and become part of the planet
                            NUM_PARTICLES_WITH_PERIAPSES_WITHIN_RADIAL_DISTANCE += 1
                            p.label = "PLANET"
                            p.assigned_body = planet_label
                            NEW_MASS_PROTOPLANET += p.mass
                            NEW_Z_ANGULAR_MOMENTUM_PROTOPLANET += p.angular_momentum_vector[2]
                        else:  # the particle is part of the disk
                            NUM_PARTICLES_IN_DISK += 1
                            p.label = "DISK"
                            p.assigned_body = None
                            NEW_MASS_DISK += p.mass
                            NEW_Z_ANGULAR_MOMENTUM_DISK += p.angular_momentum_vector[2]
                    else:  # parabolic orbit, will escape the disk
                        NUM_PARTICLES_ESCAPING += 1
                        p.label = "ESCAPE"
                        p.assigned_body = None
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
            total_angular_momentum = sum([i.angular_momentum for i in particles])
            print(
                "ITERATION: {}\n"
                "ERROR: {}\n"
                "NEW A: {}\n"
                "NUM_PARTICLES_WITHIN_RADIAL_DISTANCE: {}\n"
                "NUM_PARTICLES_WITH_PERIAPSES_WITHIN_RADIAL_DISTANCE: {}\n"
                "NUM_PARTICLES_IN_DISK: {}\n"
                "NUM_PARTICLES_ESCAPING: {}".format(iteration, error, self.a,
                                                    NUM_PARTICLES_WITHIN_RADIAL_DISTANCE,
                                                    NUM_PARTICLES_WITH_PERIAPSES_WITHIN_RADIAL_DISTANCE,
                                                    NUM_PARTICLES_IN_DISK, NUM_PARTICLES_ESCAPING)
            )
            print(
                "PROTOPLANET MASS: {} M_E ({} KG)\n"
                "DISK MASS: {} M_L ({} KG)\n"
                "ESCAPING MASS: {} M_L ({} KG)\n".format(
                    NEW_MASS_PROTOPLANET / EARTH_MASS, NEW_MASS_PROTOPLANET,
                    NEW_MASS_DISK / LUNAR_MASS, NEW_MASS_DISK,
                    NEW_MASS_ESCAPED / LUNAR_MASS, NEW_MASS_ESCAPED
                )
            )
            print(
                "TOTAL ANGULAR MOMENTUM: {} L_EM ({})\n\n".format(total_angular_momentum / L_EM, total_angular_momentum)
            )
            if self.__relative_velocity:
                self.target_velocity = [
                    statistics.mean(
                        [p.velocity_vector[0] for p in particles if p.label == "PLANET" and p.particle_name == 1]),
                    statistics.mean(
                        [p.velocity_vector[1] for p in particles if p.label == "PLANET" and p.particle_name == 1]),
                    statistics.mean(
                        [p.velocity_vector[2] for p in particles if p.label == "PLANET" and p.particle_name == 1])
                ]
            for p in particles:
                try:
                    if self.__relative_velocity:
                        p.relative_velocity_vector = [
                            p.velocity_vector[0] - self.target_velocity[0],
                            p.velocity_vector[1] - self.target_velocity[1],
                            p.velocity_vector[2] - self.target_velocity[2]
                        ]
                    p.recalculate_elements(mass_grav_body=self.mass_protoearth)
                except:
                    particles.remove(p)
            # plots.plot_eccentricity_elements(
            #     particles=particles,
            #     a=self.a,
            #     b=self.b
            # )
            # plots.plot_energies(
            #     particles=particles,
            #     a=self.a,
            #     b=self.b
            # )
            # plt.show()

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
            for p in target_removed_particles:
                p.position_vector[0] += self.earth_center[0]
                p.position_vector[1] += self.earth_center[1]
                p.position_vector[2] += self.earth_center[2]
            print("Solving impactor...")
            self.centering_resolution = 1e2
            self.centering_delta = 1.8e3
            self.earth_center = find_center(
                x=[p.position_vector[0] for p in target_removed_particles],
                y=[p.position_vector[1] for p in target_removed_particles],
                z=[p.position_vector[2] for p in target_removed_particles],
                mass=[p.mass for p in target_removed_particles],
                resolution=self.centering_resolution,
                delta_x=self.centering_delta,
                delta_y=self.centering_delta,
                delta_z=self.centering_delta,
                particle_ids=[p.particle_id for p in target_removed_particles],
                target_iron_centering=self.__center_on_target_iron
            )
            for p in target_removed_particles:
                p.position_vector[0] -= self.earth_center[0]
                p.position_vector[1] -= self.earth_center[1]
                p.position_vector[2] -= self.earth_center[2]
            impactor = self.__solve(particles=target_removed_particles, planet_label=impactor_label)
            print("Finished solving impactor!")

        if self.plot:
            plots.scatter_particles(
                x=[p.position_vector[0] for p in particles],
                y=[p.position_vector[1] for p in particles],
                tags=[p.particle_id for p in particles],
                x_label="x",
                y_label="y",
                a=self.a,
                b=self.b,
                center_plot=self.__center_plot
            )
            plots.colorcode_orbits(
                particles=particles,
                a=self.a,
                b=self.b,
                center_plot=self.__center_plot
            )
            plots.plot_eccentricities(
                particles=particles,
                a=self.a,
                b=self.b
            )

        return particles


class ParticleMapFromFiles:

    def __init__(self, path):
        self.path = path

    def read(self, time):
        particles = []
        df = pd.read_csv(self.path + "/{}.csv".format(time))
        for row in df.index:
            position_vector = [df["x"][row], df["y"][row], df["z"][row]]
            absolute_velocity_vector = [df["v_x_absolute"][row], df["v_y_absolute"][row], df["v_z_absolute"][row]]
            relative_velocity_vector = [df["v_x_relative"][row], df["v_y_relative"][row], df["v_z_relative"][row]]
            p = Particle(
                particle_name=df["particle_id"][row],
                particle_id=df["tag"][row],
                position_vector=position_vector,
                velocity_vector=absolute_velocity_vector,
                relative_velocity_vector=relative_velocity_vector,
                mass=df["mass"][row],
                density=df["density"][row],
                internal_energy=df["internal_energy"][row],
                entropy=df["entropy"][row],
                temperature=df["temperature"][row],
                mass_grav_body=df["mass_grav_body"][row],
                find_orbital_elements=False,
            )
            p.label = df["label"][row]
            p.pressure = df['pressure'][row]
            p.angular_momentum = df['angular_momentum'][row]
            p.momentum_vector = [df['total_momentum_x'][row], df['total_momentum_y'][row], df['total_momentum_z'][row]]
            particles.append(p)
        return particles
