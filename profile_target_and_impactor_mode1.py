import os
import shutil
import matplotlib.pyplot as plt
import numpy as np
from src.identify import ParticleMapFromFiles

start_time = 0
end_time = 10
interval = 1
path = "/scratch/shull4/GI2_outfiles"
output_path = "mode1_profile_plot_outputs"

if output_path in os.listdir(os.getcwd()):
    shutil.rmtree(output_path)
os.mkdir(output_path)

for time in np.arange(start_time, end_time + interval, interval):
    particle_map = ParticleMapFromFiles(path=path).read(time=time)
    target_silicate_particles = [i for i in particle_map if i.particle_id == 0]
    target_iron_particles = [i for i in particle_map if i.particle_id == 1]
    impactor_silicate_particles = [i for i in particle_map if i.particle_id == 2]
    impactor_iron_particles = [i for i in particle_map if i.particle_id == 3]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(
        [i.distance / 1000.0 for i in target_silicate_particles],
        [i.entropy for i in target_silicate_particles],
        marker="+",
        color="red",
        label="Target Silicate"
    )
    ax.scatter(
        [i.distance / 1000.0 for i in target_iron_particles],
        [i.entropy for i in target_iron_particles],
        marker="+",
        color="blue",
        label="Target Iron"
    )
    ax.scatter(
        [i.distance / 1000.0 for i in impactor_silicate_particles],
        [i.entropy for i in impactor_silicate_particles],
        marker="+",
        color="green",
        label="Impactor Silicate"
    )
    ax.scatter(
        [i.distance / 1000.0 for i in impactor_iron_particles],
        [i.entropy for i in impactor_iron_particles],
        marker="+",
        color="purple",
        label="Impactor Iron"
    )
    ax.set_xlabel("Distance from Target Center (km)")
    ax.set_ylabel("Entropy")
    ax.set_title("Time Iteration: {}".format(time))
    ax.grid()
    ax.legend(loc='lower right')
    plt.savefig(output_path + "/{}.png".format(time))
    plt.close()
