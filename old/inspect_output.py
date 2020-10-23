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
        y_values = list(self.combined_output[y_index])
        target_silicate_x = [x_values[index] for index, i in enumerate(self.ids) if (i == 0)]
        target_silicate_y = [y_values[index] for index, i in enumerate(self.ids) if i == 0]
        target_iron_x = [x_values[index] for index, i in enumerate(self.ids) if i == 1]
        target_iron_y = [y_values[index] for index, i in enumerate(self.ids) if i == 1]
        impactor_silicate_x = [x_values[index] for index, i in enumerate(self.ids) if (i == 2)]
        impactor_silicate_y = [y_values[index] for index, i in enumerate(self.ids) if i == 2]
        impactor_iron_x = [x_values[index] for index, i in enumerate(self.ids) if i == 3]
        impactor_iron_y = [y_values[index] for index, i in enumerate(self.ids) if i == 3]
        ax = plt.figure().add_subplot(111)
        ax.scatter(target_silicate_x, target_silicate_y, marker="+", color="red", label="target silicate")
        ax.scatter(target_iron_x, target_iron_y, marker="+", color="blue", label="target iron")
        ax.scatter(impactor_silicate_x, impactor_silicate_y, marker="+", color="yellow", label="impactor silicate")
        ax.scatter(impactor_iron_x, impactor_iron_y, marker="+", color="green", label="impactor iron")
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.legend()
        ax.grid()
        plt.show()


m = Inspect(num_processes=20, time=6, output_path="/Users/scotthull/Desktop/GI")
m.plot(9, 11, "Density", "Pressure")


