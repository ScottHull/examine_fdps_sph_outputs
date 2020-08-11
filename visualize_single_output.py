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

    def plot(self, dimension):
        x, y, z, mass, particle_id = self.__return_data()
        if dimension == 2:
            ax = plt.figure().add_subplot(111)
            ax.scatter(x, y, color='black', marker="+")
        else:
            ax = Axes3D(plt.figure())
            ax.scatter(x, y, z, color="black", marker="+")
            ax.set_zlabel("z")
        ax.set_xlabel("x")
        ax.set_xlabel("y")

        return ax


v = Visualize(output_path="/Users/scotthull/Desktop/imp.dat")
ax = v.plot(dimension=2)
plt.show()

