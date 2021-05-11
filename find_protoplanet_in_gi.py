import os
from math import sqrt
import pandas as pd
import numpy as np
from src.identify import ParticleMap
from src.combine import CombineFile
from src.report import make_report


def calc_escape_velocity(mass, radius):
    G = 6.674 * 10 ** -11
    return sqrt((2 * G * mass) / radius)


time = 2000
number_processes = 100
path_to_outputs = "/scratch/shull4/gi"

combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
f = os.getcwd() + "/merged_{}.dat".format(time)
pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
particle_map = pm.solve()

mass = sum([p.mass for p in particle_map if p.label == "PLANET"])
escape_velocity = calc_escape_velocity(mass=mass, radius=pm.a)

print(
    "Mass: {}\nRadius: {}\nEscape Velocity: {}".format(mass, pm.a, escape_velocity)
)
