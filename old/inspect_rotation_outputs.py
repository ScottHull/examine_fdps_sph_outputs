import os
from math import sqrt
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Visualize:

    def __init__(self, output_path, num_files, start_from=0):
        self.output_path = output_path
        self.num_files = num_files
        self.start_file = start_from
        self.curr_file = start_from
        self.file_format = "results.{}_00001_00000.dat"

    def __get_filename(self):
        return self.output_path + "/" + self.file_format.format(str(self.curr_file).zfill(5))

    def __read_sph_file(self):
        df = pd.read_csv(self.__get_filename(), sep='\t', skiprows=2, header=None)
        x = df[3]
        y = df[4]
        z = df[5]
        v_x = df[6]
        v_y = df[7]
        v_z = df[8]
        return x, y, z, v_x, v_y, v_z

    def __make_scene(self, ax):
        x, y, z, v_x, v_y, v_z = self.__read_sph_file()
        r = [sqrt(i ** 2 + j ** 2) for i, j in zip(x, y)]
        v = [sqrt(i ** 2 + j ** 2) for i, j in zip(v_x, v_y)]
        v_omega = [1E-4 * i for i in r]
        ax.scatter(r, v, c='black')
        ax.scatter(r, v_omega, c='red')

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i in range(self.start_file, self.num_files + 1):
            self.__make_scene(ax=ax)
        ax.set_xlabel('sqrt(x^2 + y^2)')
        ax.set_ylabel('sqrt(v_x^2 + v_y^2)')
        plt.show()


v = Visualize(output_path="/Users/scotthull/Documents/FDPS_SPH/test2", num_files=1000, start_from=950)
v.plot()
