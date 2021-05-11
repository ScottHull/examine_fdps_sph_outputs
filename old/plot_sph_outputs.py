import os
from math import sqrt
from copy import copy
import shutil
import pandas as pd
import moviepy.editor as mpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Output:

    def __init__(self, path_to_output):
        self.df = pd.read_csv(path_to_output, sep='\t', skiprows=2, header=None)
        self.id = self.df[0]
        self.tag = self.df[1]
        self.x = self.df[3]
        self.y = self.df[4]
        self.z = self.df[5]
        self.radius = [sqrt(x ** 2 + y ** 2 + z ** 2) / 1000 for x, y, z in zip(self.x, self.y, self.z)]
        self.target_silicate_radius = [self.radius[index] for index, i in enumerate(self.tag) if i == 0]
        self.target_iron_radius = [self.radius[index] for index, i in enumerate(self.tag) if i == 1]
        self.impactor_silicate_radius = [self.radius[index] for index, i in enumerate(self.tag) if i == 2]
        self.impactor_iron_radius = [self.radius[index] for index, i in enumerate(self.tag) if i == 3]
        self.pressure = self.df[11]
        self.density = self.df[9]
        self.energy = self.df[10]
        # self.entropy = self.df[13]
        # self.soundspeed = self.df[15]
        # self.dt = self.df[16]

        for index, i in enumerate(self.id):
            p = self.pressure[index] * 10 ** -9
            d = self.density[index]
            u = self.energy[index]
            if p > 1700:
                print(i, p, d, u)

    def plot_pressure(self):
        target_silicate_pressure = [self.pressure[index] * 10 ** -9 for index, i in enumerate(self.tag) if i == 0]
        target_iron_pressure = [self.pressure[index] * 10 ** -9 for index, i in enumerate(self.tag) if i == 1]
        impactor_silicate_pressure = [self.pressure[index] * 10 ** -9 for index, i in enumerate(self.tag) if i == 2]
        impactor_iron_pressure = [self.pressure[index] * 10 ** -9 for index, i in enumerate(self.tag) if i == 3]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(self.target_silicate_radius, target_silicate_pressure, color='blue', marker="+",
                   label="Target Silicate")
        ax.scatter(self.target_iron_radius, target_iron_pressure, color='red', marker="+", label="Target Iron")
        ax.scatter(self.impactor_silicate_radius, impactor_silicate_pressure, color='green', marker="+",
                   label="Impactor Silicate")
        ax.scatter(self.impactor_iron_radius, impactor_iron_pressure, color='purple', marker="+", label="Impactor Iron")
        ax.set_xlabel("Radius (km)")
        ax.set_ylabel("Pressure (GPa)")
        ax.set_title("Pressure")
        ax.legend(loc='upper left')
        ax.grid()
        plt.show()

    def plot_density(self):
        target_silicate_density = [self.density[index] for index, i in enumerate(self.tag) if i == 0]
        target_iron_density = [self.density[index] for index, i in enumerate(self.tag) if i == 1]
        impactor_silicate_density = [self.density[index] for index, i in enumerate(self.tag) if i == 2]
        impactor_iron_density = [self.density[index] for index, i in enumerate(self.tag) if i == 3]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(self.target_silicate_radius, target_silicate_density, color='blue', marker="+",
                   label="Target Silicate")
        ax.scatter(self.target_iron_radius, target_iron_density, color='red', marker="+", label="Target Iron")
        ax.scatter(self.impactor_silicate_radius, impactor_silicate_density, color='green', marker="+",
                   label="Impactor Silicate")
        ax.scatter(self.impactor_iron_radius, impactor_iron_density, color='purple', marker="+", label="Impactor Iron")
        ax.set_xlabel("Radius (km)")
        ax.set_ylabel("Density (kg/m3)")
        ax.set_title("Density")
        ax.legend(loc='upper left')
        ax.grid()
        plt.show()

    def plot_energy(self):
        target_silicate_energy = [self.energy[index] for index, i in enumerate(self.tag) if i == 0]
        target_iron_energy = [self.energy[index] for index, i in enumerate(self.tag) if i == 1]
        impactor_silicate_energy = [self.energy[index] for index, i in enumerate(self.tag) if i == 2]
        impactor_iron_energy = [self.energy[index] for index, i in enumerate(self.tag) if i == 3]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(self.target_silicate_radius, target_silicate_energy, color='blue', marker="+",
                   label="Target Silicate")
        ax.scatter(self.target_iron_radius, target_iron_energy, color='red', marker="+", label="Target Iron")
        ax.scatter(self.impactor_silicate_radius, impactor_silicate_energy, color='green', marker="+",
                   label="Impactor Silicate")
        ax.scatter(self.impactor_iron_radius, impactor_iron_energy, color='purple', marker="+", label="Impactor Iron")
        ax.set_xlabel("Radius (km)")
        ax.set_ylabel("Energy (J/kg)")
        ax.set_title("Energy")
        ax.legend(loc='upper left')
        ax.grid()
        plt.show()

    # def plot_soundspeed(self):
    #     target_silicate_soundspeed = [self.soundspeed[index] for index, i in enumerate(self.tag) if i == 0]
    #     target_iron_soundspeed = [self.soundspeed[index] for index, i in enumerate(self.tag) if i == 1]
    #     impactor_silicate_soundspeed = [self.soundspeed[index] for index, i in enumerate(self.tag) if i == 2]
    #     impactor_iron_soundspeed = [self.soundspeed[index] for index, i in enumerate(self.tag) if i == 3]
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111)
    #     ax.scatter(self.target_silicate_radius, target_silicate_soundspeed, color='blue', marker="+",
    #                label="Target Silicate")
    #     ax.scatter(self.target_iron_radius, target_iron_soundspeed, color='red', marker="+", label="Target Iron")
    #     ax.scatter(self.impactor_silicate_radius, impactor_silicate_soundspeed, color='green', marker="+",
    #                label="Impactor Silicate")
    #     ax.scatter(self.impactor_iron_radius, impactor_iron_soundspeed, color='purple', marker="+", label="Impactor Iron")
    #     ax.set_xlabel("Radius (km)")
    #     ax.set_ylabel("Soundspeed")
    #     ax.set_title("Soundspeed")
    #     ax.legend(loc='upper left')
    #     ax.grid()
    #     plt.show()
    #
    # def plot_dt(self):
    #     target_silicate_dt = [self.dt[index] for index, i in enumerate(self.tag) if i == 0]
    #     target_iron_dt = [self.dt[index] for index, i in enumerate(self.tag) if i == 1]
    #     impactor_silicate_dt = [self.dt[index] for index, i in enumerate(self.tag) if i == 2]
    #     impactor_iron_dt = [self.dt[index] for index, i in enumerate(self.tag) if i == 3]
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111)
    #     ax.scatter(self.target_silicate_radius, target_silicate_dt, color='blue', marker="+",
    #                label="Target Silicate")
    #     ax.scatter(self.target_iron_radius, target_iron_dt, color='red', marker="+", label="Target Iron")
    #     ax.scatter(self.impactor_silicate_radius, impactor_silicate_dt, color='green', marker="+",
    #                label="Impactor Silicate")
    #     ax.scatter(self.impactor_iron_radius, impactor_iron_dt, color='purple', marker="+", label="Impactor Iron")
    #     ax.set_xlabel("Radius (km)")
    #     ax.set_ylabel("dt")
    #     ax.set_title("Timestep (dt)")
    #     ax.legend(loc='upper left')
    #     ax.grid()
    #     plt.show()


m = Output(path_to_output="/Users/scotthull/Desktop/imp_new_gcc.dat")
m.plot_density()
m.plot_pressure()
m.plot_energy()
