import os
import shutil
from math import sqrt
import numpy as np
import pandas as pd
import moviepy.editor as mpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.combine import CombineFile

start_time = 0
end_time = 1000
interval = 1
number_processes = 100
path = "/scratch/shull4/impactor"
generic_output_path = "/scratch/shull4/single_body_animate"
internal_energy_path = "/scratch/shull4/single_body_animate_internal_energy"
paths = [generic_output_path, internal_energy_path]
fnames = ["single_body_output.mp4", "single_body_output_internal_energy.mp4"]

for output_path in paths:
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)


def animate(start_time, end_time, interval, output_path, file_name="single_body_output.mp4"):
    frames = [output_path + "/{}.png".format(i) for i in np.arange(start_time, end_time + interval, interval)]
    animation = mpy.ImageSequenceClip(frames, fps=30, load_images=True)
    animation.write_videofile(file_name, fps=30)


for time in np.arange(start_time, end_time + interval, interval):
    combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path).combine()
    sph_file = os.getcwd() + "/merged_{}.dat".format(time)
    df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")
    os.remove(sph_file)
    tag, x, y, z, energy, entropy = df[1], df[3], df[4], df[5], df[10], df[13]

    silicate_x = [i for index, i in enumerate(x) if tag[index] == 0]
    silicate_y = [i for index, i in enumerate(y) if tag[index] == 0]
    silicate_z = [i for index, i in enumerate(z) if tag[index] == 0]
    silicate_internal_energy = [i for index, i in enumerate(energy) if tag[index] == 0]
    silicate_entropy = [i for index, i in enumerate(entropy) if tag[index] == 0]
    silicate_distances = [sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2) / 1000.0 for i in
                          zip(silicate_x, silicate_y, silicate_z)]
    iron_x = [i for index, i in enumerate(x) if tag[index] == 1]
    iron_y = [i for index, i in enumerate(y) if tag[index] == 1]
    iron_z = [i for index, i in enumerate(z) if tag[index] == 1]
    iron_internal_energy = [i for index, i in enumerate(energy) if tag[index] == 1]
    iron_entropy = [i for index, i in enumerate(entropy) if tag[index] == 1]
    iron_distances = [sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2) / 1000.0 for i in zip(iron_x, iron_y, iron_z)]

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    ax_silicate = fig.add_subplot(121, projection='3d')
    ax_iron = fig.add_subplot(122, projection='3d')
    ax_silicate.scatter(
        silicate_x,
        silicate_y,
        silicate_z,
        color="red",
        marker="+"
    )
    ax_iron.scatter(
        iron_x,
        iron_y,
        iron_z,
        color="blue",
        marker="+"
    )
    ax_silicate.set_title("Silicate")
    ax_iron.set_title("Iron")
    ax_silicate.set_xlim(-8000e3, 8000e3)
    ax_silicate.set_ylim(-8000e3, 8000e3)
    ax_silicate.set_zlim(-8000e3, 8000e3)
    ax_iron.set_xlim(-8000e3, 8000e3)
    ax_iron.set_ylim(-8000e3, 8000e3)
    ax_iron.set_zlim(-8000e3, 8000e3)
    plt.tight_layout()
    plt.savefig(generic_output_path + "/{}.png".format(time), format="png")
    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    sc = ax.scatter(
        silicate_distances + iron_distances,
        silicate_internal_energy + iron_internal_energy,
        c=silicate_entropy + iron_entropy,
        marker="+"
    )
    cbar = plt.colorbar(sc)
    cbar.set_label("Entropy")
    ax.set_xlabel("Distance from Target Center")
    ax.set_ylabel("Internal Energy")
    ax.set_title("Iteration: {}".format(time))
    ax.grid()

for index, output_path in enumerate(paths):
    fname = fnames[index]
    animate(start_time=start_time, end_time=end_time, interval=interval, output_path=output_path, file_name=fname)
    shutil.rmtree(output_path)
