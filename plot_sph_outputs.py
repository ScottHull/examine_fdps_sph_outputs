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
        self.id = self.df[1]
        self.x = self.df[3]
        self.y = self.df[4]
        self.z = self.df[5]
        self.radius = [sqrt(x**2 + y**2 + z**2) / 1000 for x, y, z in zip(self.x, self.y, self.z)]
        self.pressure = self.df[11]
        self.density = self.df[9]
        self.energy = self.df[10]
        self.entropy = self.df[13]

    def plot_pressure(self):
        silicate_pressure = [self.pressure[index] * 10**-9 for index, i in enumerate(self.id) if i == 0]
        iron_pressure = [self.pressure[index] * 10**-9 for index, i in enumerate(self.id) if i == 1]
        silicate_radius = [self.radius[index] for index, i in enumerate(self.id) if i == 0]
        iron_radius = [self.radius[index] for index, i in enumerate(self.id) if i == 1]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(silicate_radius, silicate_pressure, color='blue', marker="+", label="silicate")
        ax.scatter(iron_radius, iron_pressure, color='red', marker="+", label="iron")
        ax.set_xlabel("Radius (km)")
        ax.set_ylabel("Pressure (GPa)")
        ax.set_title("Pressure")
        ax.legend(loc='upper left')
        ax.grid()
        plt.show()
        
    def plot_density(self):
        silicate_density = [self.density[index] for index, i in enumerate(self.id) if i == 0]
        iron_density = [self.density[index] for index, i in enumerate(self.id) if i == 1]
        silicate_radius = [self.radius[index] for index, i in enumerate(self.id) if i == 0]
        iron_radius = [self.radius[index] for index, i in enumerate(self.id) if i == 1]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(silicate_radius, silicate_density, color='blue', marker="+", label="silicate")
        ax.scatter(iron_radius, iron_density, color='red', marker="+", label="iron")
        ax.set_xlabel("Radius (km)")
        ax.set_ylabel("Density (kg/m3)")
        ax.set_title("Density")
        ax.legend(loc='upper left')
        ax.grid()
        plt.show()
        
    def plot_energy(self):
        silicate_energy = [self.energy[index] for index, i in enumerate(self.id) if i == 0]
        iron_energy = [self.energy[index] for index, i in enumerate(self.id) if i == 1]
        silicate_radius = [self.radius[index] for index, i in enumerate(self.id) if i == 0]
        iron_radius = [self.radius[index] for index, i in enumerate(self.id) if i == 1]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(silicate_radius, silicate_energy, color='blue', marker="+", label="silicate")
        ax.scatter(iron_radius, iron_energy, color='red', marker="+", label="iron")
        ax.set_xlabel("Radius (km)")
        ax.set_ylabel("Energy (J/kg)")
        ax.set_title("Energy")
        ax.legend(loc='upper left')
        ax.grid()
        plt.show()

m = Output(path_to_output="/Users/scotthull/Documents/FDPS_SPH/test2/results.01000_00001_00000.dat")
m.plot_density()
m.plot_pressure()
m.plot_energy()
