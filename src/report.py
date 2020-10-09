import os
import pandas as pd


def make_report(particles, time, to_directory="/particle_outputs"):
    if not os.path.exists(os.getcwd() + to_directory):
        os.mkdir(os.getcwd() + to_directory)
    path = os.getcwd() + to_directory
    pd.DataFrame({
        "particle_id": [p.particle_name for p in particles],
        "tag": [p.particle_id for p in particles],
        "x": [p.position_vector[0] for p in particles],
        "y": [p.position_vector[1] for p in particles],
        "z": [p.position_vector[2] for p in particles],
        "v_x": [p.velocity_vector[0] for p in particles],
        "v_y": [p.velocity_vector[1] for p in particles],
        "v_z": [p.velocity_vector[2] for p in particles],
        "eccentricity": [p.eccentricity for p in particles],
        "orbital_energy": [p.orbital_energy for p in particles]
    }).to_csv(path + "/{}.csv".format(time))
