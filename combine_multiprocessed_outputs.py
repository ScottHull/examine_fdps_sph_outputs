import pandas as pd

class CombineFile:

    def __init__(self, num_processes, time, output_path):
        self.num_processes = num_processes
        self.time = time
        self.output_path = output_path
        self.file_format = "results.{}_0000{}_0000{}.dat"
        self.curr_process = 0

    def __get_filename(self):
        return self.output_path + "/" + self.file_format.format(str(self.time).zfill(5), self.num_processes,
                                                                self.curr_process)

    def __read_sph_file(self):
        df = pd.read_csv(self.__get_filename(), sep='\t', skiprows=2, header=None)
        return df

    def combine(self):
        dfs = []
        for proc in range(0, self.num_processes, 1):
            self.curr_process = proc
            dfs.append(self.__read_sph_file())
        merged_df = pd.concat(dfs)
        merged_df.to_csv("merged_{}.dat".format(self.time))

