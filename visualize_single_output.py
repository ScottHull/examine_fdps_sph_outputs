import os
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Visualize:

    def __init__(self, output_path):
        self.output_path = output_path

    def __return_data(self):
        df = pd.read_csv(self.output_path, sep='\t', skiprows=2, header=None)
        particle_id = df[1]
        mass = df[2]
        x = df[3]
        y = df[4]
        z = df[5]

        return x, y, z, mass, particle_id

    def __return_color_list(self, header_index):
        df = pd.read_csv(self.output_path, sep='\t', skiprows=2, header=None)
        return df[header_index]

    def plot(self, dimension, color=None, color_label=""):
        x, y, z, mass, particle_id = self.__return_data()
        if dimension == 2:
            ax = plt.figure().add_subplot(111)
            if color is None:
                ax.scatter(x, y, color='black', marker="+")
            else:
                sc = ax.scatter(x, y, c=self.__return_color_list(header_index=color), marker="+")
                cbar = plt.colorbar(sc, ax=ax)
                cbar.set_label(color_label)
        else:
            ax = Axes3D(plt.figure())
            ax.scatter(x, y, z, color="black", marker="+")
            ax.set_zlabel("z")
        ax.set_xlabel("x")
        ax.set_xlabel("y")

        return ax


v = Visualize(output_path="/Users/scotthull/Desktop/merged_212.dat")
ax = v.plot(dimension=2, color=15, color_label="Soundspeed")
plt.show()

