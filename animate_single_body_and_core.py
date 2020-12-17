import os
import shutil
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
output_path = "/scratch/shull4/single_body_animate"

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
    tag, x, y, z = df[1], df[3], df[4], df[5]

    silicate_x = [i for index, i in enumerate(x) if tag[index] == 0]
    silicate_y = [i for index, i in enumerate(y) if tag[index] == 0]
    silicate_z = [i for index, i in enumerate(z) if tag[index] == 0]
    iron_x = [i for index, i in enumerate(x) if tag[index] == 1]
    iron_y = [i for index, i in enumerate(y) if tag[index] == 1]
    iron_z = [i for index, i in enumerate(z) if tag[index] == 1]

    fig = plt.figure()
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
    plt.savefig(output_path + "/{}.png".format(time), format="png")
    plt.close()

animate(start_time=start_time, end_time=end_time, interval=interval, output_path=output_path)
shutil.rmtree(output_path)
