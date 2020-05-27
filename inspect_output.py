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

    def plot(self, x_index, y_index, x_label, y_label):
        ax = plt.figure().add_subplot(111)
        ax.scatter(self.combined_output[x_index], self.combined_output[y_index], marker="+", color="black")
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.grid()
        plt.show()


m = Inspect(num_processes=1, time=6, output_path="/Users/scotthull/Desktop/GI")
m.plot(9, 11, "Density", "Pressure")


