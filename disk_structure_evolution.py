import os
import pandas as pd
import numpy as np
from src.identify import ParticleMap
from src.combine import CombineFile
from src.report import make_report
import src.plots as plots
from src.structure import Structure
import matplotlib.pyplot as plt
import shutil

start_time = 50
end_time = 1000
interval = 50
number_processes = 100

path_to_outputs = "/scratch/shull4/GI2/"
eccentricity_plot_path = os.getcwd() + "/eccentricity"
disk_structure_path = os.getcwd() + "/structure"
disk_structure_eccentricity_path = os.getcwd() + "/eccentricity_structure"
vmf_path = os.getcwd() + "/vmf"

paths = [eccentricity_plot_path, disk_structure_path, disk_structure_eccentricity_path]
for i in paths:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

for time in np.arange(start_time, end_time + interval, interval):
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
    particle_map = pm.solve()
    s = Structure(disk_particles=particle_map, phase="duniteS2")
    # vmf = s.calc_vapor_mass_fraction()
    # make_report(particles=particle_map, time=time)
    os.remove(f)
    fig = plots.plot_eccentricities(particles=particle_map, a=pm.a, b=pm.b)
    fig.savefig(eccentricity_plot_path + "/{}.png".format(time), format='png')
    plt.close()
    fig = plots.scatter_particles(
                x=[p.position_vector[0] for p in particle_map],
                y=[p.position_vector[1] for p in particle_map],
                tags=[p.particle_id for p in particle_map],
                x_label="x",
                y_label="y",
                a=pm.a,
                b=pm.b,
                center_plot=True
            )
    fig.savefig(disk_structure_path + "/{}.png".format(time), format='png')
    plt.close()
    fig = plots.colorcode_orbits(
                particles=particle_map,
                a=pm.a,
                b=pm.b,
                center_plot=True
            )
    fig.savefig(disk_structure_eccentricity_path + "/{}.png".format(time), format='png')
    plt.close()
    fig = plots.plot_vfm(
        phase_curve_1_x=s.phase_df['entropy_sol_liq'],
        phase_curve_1_y=s.phase_df['temperature'],
        phase_curve_2_x=s.phase_df['entropy_vap'],
        phase_curve_2_y=s.phase_df['temperature'],
        particles_x=[p.entropy for p in particle_map],
        particles_y=[p.temperature for p in particle_map],
        particles_color=[p.distance for p in particle_map],
        xlabel="Entropy (S)",
        ylabel="Temperature (T [deg K])",
        cbar_label="Radial Distance From Target Center",
        phase_curve_1_label="sol-liq",
        phase_curve_2_label="liq-gas"
    )
    fig.savefig(vmf_path + "/{}.png".format(time), format='png')
    plt.close()

plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=eccentricity_plot_path,
              filename="eccentricities.mp4", fps=5)
plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=disk_structure_path,
              filename="structure.mp4", fps=5)
plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=disk_structure_eccentricity_path,
              filename="structure_eccentricities.mp4", fps=5)
plots.animate(start_time=start_time, end_time=end_time, interval=interval, path=vmf_path,
              filename="vmf.mp4", fps=5)

