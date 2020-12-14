import pandas as pd
import numpy as np

# this code works reports.py outputs for pre-profiled simulations
start_time = 0
end_time = 5000
interval = 1
path = "/scratch/shull4/GI_outfiles"

for time in np.arange(start_time, end_time + interval, interval):
    particles = pd.read_csv(path + "/{}.csv".format(time))
