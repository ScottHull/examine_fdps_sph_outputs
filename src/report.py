import os
import pandas as pd


def make_report(particles, time, to_directory="/particle_outputs"):
    pd.DataFrame({
        "particle_id": [p.particle_name for p in particles],
        "tag": [p.particle_id for p in particles],
        "label": [p.label for p in particles],
        "x": [p.position_vector[0] for p in particles],
        "y": [p.position_vector[1] for p in particles],
        "z": [p.position_vector[2] for p in particles],
        "v_x_absolute": [p.velocity_vector[0] for p in particles],
        "v_y_absolute": [p.velocity_vector[1] for p in particles],
        "v_z_absolute": [p.velocity_vector[2] for p in particles],
        "v_x_relative": [p.relative_velocity_vector[0] for p in particles],
        "v_y_relative": [p.relative_velocity_vector[1] for p in particles],
        "v_z_relative": [p.relative_velocity_vector[2] for p in particles],
        "eccentricity": [p.eccentricity for p in particles],
        "orbital_energy": [p.orbital_energy for p in particles],
        "mass": [p.mass for p in particles],
        "internal_energy": [p.internal_energy for p in particles],
        "density": [p.density for p in particles],
        "entropy": [p.entropy for p in particles],
        "temperature": [p.temperature for p in particles],
        "mass_grav_body": [p.mass_grav_body for p in particles],
        "pressure": [p.pressure for p in particles],
        "angular_momentum": [p.angular_momentum for p in particles]
    }).to_csv(to_directory + "/{}.csv".format(time))
