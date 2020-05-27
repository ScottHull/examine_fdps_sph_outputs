import pandas as pd
import csv
import os
import matplotlib.pyplot as plt
from combine_multiprocessed_outputs import CombineFile


class Inspect:

    def __init__(self, num_processes, time, output_path):
        self.num_processes = num_processes
        self.time = time
        self.output_path = output_path
        self.file_format = "results.{}_{}_{}.dat"
        self.curr_process = 0
        self.combined_output = CombineFile(num_processes=self.num_processes, time=self.time,
                                           output_path=self.output_path).combine_df()
        self.ids = list(self.combined_output[1])

    def plot(self, x_index, y_index, x_label, y_label):
        x_values = list(self.combined_output[x_index])
        silicate_x = [self.combined_output[x_index][index] for index, i in enumerate(self.ids) if i == 0]
        silicate_y = [self.combined_output[y_index][index] for index, i in enumerate(self.ids) if i == 0]
        iron_x = [self.combined_output[x_index][index] for index, i in enumerate(self.ids) if i == 1]
        iron_y = [self.combined_output[y_index][index] for index, i in enumerate(self.ids) if i == 1]
        ax = plt.figure().add_subplot(111)
        ax.scatter(silicate_x, silicate_y, marker="+", color="red", label="silicate")
        ax.scatter(iron_x, iron_y, marker="+", color="blue", label="iron")
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.legend()
        ax.grid()
        plt.show()


m = Inspect(num_processes=20, time=12, output_path="/Users/scotthull/Desktop/GI")
m.plot(9, 11, "Density", "Pressure")


