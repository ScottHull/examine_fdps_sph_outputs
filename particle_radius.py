from math import sqrt
import pandas as pd
import matplotlib.pyplot as plt



class RadiusPlot:

    def __init__(self, path, focus_timestep, num_processes=1, focus_process=None):
        self.path = path
        self.num_processes = num_processes
        self.focus_process = focus_process
        self.curr_process = 0
        self.file_format = "results.{}_0000{}_0000{}.dat"
        self.focus_timestep = focus_timestep


    def __get_filename(self):
        return self.path + "/" + self.file_format.format(str(self.focus_timestep).zfill(5), self.num_processes,
                                                                self.curr_process)

    def __read_sph_file(self):
        df = pd.read_csv(self.__get_filename(), sep='\t', skiprows=2, header=None)
        index_id = df[0]
        particle_id = df[1]
        x = df[3]
        y = df[4]
        z = df[5]

        return index_id, particle_id, x, y, z

    def run(self):
        colors = ['red', 'green', 'blue', 'purple', 'orange', 'black']
        fig = plt.figure()
        fig2 = plt.figure()
        ax = fig.add_subplot(111)
        ax2 = fig2.add_subplot(111)
        ax.set_xlabel("Index")
        ax.set_ylabel("Radius")
        ax.set_title("Index vs. Radius (SPH Outputs)")
        ax.grid()
        ax2.set_xlabel("Index")
        ax2.set_ylabel("Radius")
        ax2.set_title("Index vs. Radius (SPH Outputs)")
        ax2.grid()
        for proc in range(0, self.num_processes, 1):
            c = colors[proc]
            self.curr_process = proc
            index_id, p_id, x, y, z = self.__read_sph_file()
            silicate_index_id = [i for index, i in enumerate(index_id) if p_id[index] == 0]
            silicate_radius = [sqrt(i ** 2 + j ** 2 + k ** 2) for p, i, j, k in zip(p_id, x, y, z) if p == 0]
            ax.scatter(silicate_index_id, silicate_radius, marker="+", color=c, label="Process: {} (silicate)".format(self.curr_process))
            iron_index_id = [i for index, i in enumerate(index_id) if p_id[index] == 1]
            iron_radius = [sqrt(i ** 2 + j ** 2 + k ** 2) for p, i, j, k in zip(p_id, x, y, z) if p == 1]
            ax.scatter(iron_index_id, iron_radius, marker="o", color=c, label="Process: {} (iron)".format(self.curr_process))

            ax2.scatter(iron_index_id, iron_radius, marker="+", color='blue', label="Iron")
            ax2.scatter(silicate_index_id, silicate_radius, marker="+", color='red', label="Silicate")

        ax.legend(loc='lower left')
        ax2.legend(loc='lower left')
        plt.show()

path = "/Users/scotthull/Desktop/mpi_bugs/2_processes_output"
m = RadiusPlot(
    path=path,
    num_processes=2,
    focus_timestep=0
)

m.run()
