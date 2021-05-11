import os
import shutil
from math import sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import moviepy.editor as mpy
from src.identify import ParticleMapFromFiles

start_time = 0
end_time = 10
interval = 1
path = "/scratch/shull4/GI_outfiles"
output_path = "mode1_profile_plot_outputs"
positions_output_path = "mode1_positions_profile_plot_outputs"
both_output_path = "mode1_positions_both_profile_plot_outputs"
combined_path = "mode1_combined_plot_outputs"
paths = [output_path, both_output_path, positions_output_path, combined_path]

for o in paths:
    if o in os.listdir(os.getcwd()):
        shutil.rmtree(o)
    os.mkdir(o)
    os.mkdir(o + "/target")
    os.mkdir(o + "/impactor")


def animate(start_time, end_time, interval, output_path, file_name="output.mp4"):
    frames = [output_path + "/{}.png".format(i) for i in np.arange(start_time, end_time + interval, interval)]
    animation = mpy.ImageSequenceClip(frames, fps=2, load_images=True)
    animation.write_videofile(file_name, fps=2)


for time in np.arange(start_time, end_time + interval, interval):
    particle_map = ParticleMapFromFiles(path=path).read(time=time)
    target_silicate_particles = [i for i in particle_map if i.particle_id == 0]
    target_iron_particles = [i for i in particle_map if i.particle_id == 1]
    impactor_silicate_particles = [i for i in particle_map if i.particle_id == 2]
    impactor_iron_particles = [i for i in particle_map if i.particle_id == 3]
    target_silicate_distances = [
        sqrt(i.position_vector[0] ** 2 + i.position_vector[1] ** 2 + i.position_vector[2] ** 2) / 1000.0 for i in
        target_silicate_particles]
    target_iron_distances = [
        sqrt(i.position_vector[0] ** 2 + i.position_vector[1] ** 2 + i.position_vector[2] ** 2) / 1000.0 for i in
        target_iron_particles]
    impactor_silicate_distances = [
        sqrt(i.position_vector[0] ** 2 + i.position_vector[1] ** 2 + i.position_vector[2] ** 2) / 1000.0 for i in
        impactor_silicate_particles]
    impactor_iron_distances = [
        sqrt(i.position_vector[0] ** 2 + i.position_vector[1] ** 2 + i.position_vector[2] ** 2) / 1000.0 for i in
        impactor_iron_particles]

    target_fig = plt.figure()
    target_fig.set_size_inches(18.5, 10.5)
    target_silicate = target_fig.add_subplot(121)
    target_iron = target_fig.add_subplot(122)
    sc_silicate = target_silicate.scatter(
        target_silicate_distances,
        [i.internal_energy for i in target_silicate_particles],
        c=[i.entropy for i in target_silicate_particles]
    )
    cbar_silicate = plt.colorbar(sc_silicate)
    cbar_silicate.set_label("Entropy")
    target_silicate.set_xlabel("Distance from Target Center (km)")
    target_silicate.set_ylabel("Internal Energy")
    target_silicate.set_title("Silicate: Iteration: {}".format(time))
    target_silicate.grid()
    sc_iron = target_iron.scatter(
        target_iron_distances,
        [i.internal_energy for i in target_iron_particles],
        c=[i.entropy for i in target_iron_particles]
    )
    cbar_iron = plt.colorbar(sc_iron)
    cbar_iron.set_label("Entropy")
    target_iron.set_xlabel("Distance from Target Center (km)")
    target_iron.set_ylabel("Internal Energy")
    target_iron.set_title("Iron: Iteration: {}".format(time))
    target_iron.grid()
    plt.tight_layout()
    plt.savefig(output_path + "/target/{}.png".format(time), format="png")
    plt.close()

    impactor_fig = plt.figure()
    impactor_fig.set_size_inches(18.5, 10.5)
    impactor_silicate = impactor_fig.add_subplot(121)
    impactor_iron = impactor_fig.add_subplot(122)
    sc_silicate = impactor_silicate.scatter(
        impactor_silicate_distances,
        [i.internal_energy for i in impactor_silicate_particles],
        c=[i.entropy for i in impactor_silicate_particles]
    )
    cbar_silicate = plt.colorbar(sc_silicate)
    cbar_silicate.set_label("Entropy")
    impactor_silicate.set_xlabel("Distance from Target Center (km)")
    impactor_silicate.set_ylabel("Internal Energy")
    impactor_silicate.set_title("Silicate: Iteration: {}".format(time))
    impactor_silicate.grid()
    sc_iron = impactor_iron.scatter(
        impactor_iron_distances,
        [i.internal_energy for i in impactor_iron_particles],
        c=[i.entropy for i in impactor_iron_particles]
    )
    cbar_iron = plt.colorbar(sc_iron)
    cbar_iron.set_label("Entropy")
    impactor_iron.set_xlabel("Distance from impactor Center (km)")
    impactor_iron.set_ylabel("Internal Energy")
    impactor_iron.set_title("Iron: Iteration: {}".format(time))
    impactor_iron.grid()
    plt.tight_layout()
    plt.savefig(output_path + "/impactor/{}.png".format(time), format="png")
    plt.close()

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    ax_silicate = fig.add_subplot(121, projection='3d')
    ax_iron = fig.add_subplot(122, projection='3d')
    ax_silicate.scatter(
        [i.position_vector[0] for i in target_silicate_particles],
        [i.position_vector[1] for i in target_silicate_particles],
        [i.position_vector[2] for i in target_silicate_particles],
        color="red",
        marker="+"
    )
    ax_iron.scatter(
        [i.position_vector[0] for i in target_iron_particles],
        [i.position_vector[1] for i in target_iron_particles],
        [i.position_vector[2] for i in target_iron_particles],
        color="blue",
        marker="+"
    )
    ax_silicate.set_title("Target Silicate: Iteration {}".format(time))
    ax_iron.set_title("Target Iron: Iteration {}".format(time))
    ax_silicate.set_xlim(-8000e3, 8000e3)
    ax_silicate.set_ylim(-8000e3, 8000e3)
    ax_silicate.set_zlim(-8000e3, 8000e3)
    ax_iron.set_xlim(-8000e3, 8000e3)
    ax_iron.set_ylim(-8000e3, 8000e3)
    ax_iron.set_zlim(-8000e3, 8000e3)
    plt.tight_layout()
    plt.savefig(positions_output_path + "/target/{}.png".format(time), format="png")
    plt.close()

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    ax_silicate = fig.add_subplot(121, projection='3d')
    ax_iron = fig.add_subplot(122, projection='3d')
    ax_silicate.scatter(
        [i.position_vector[0] for i in impactor_silicate_particles],
        [i.position_vector[1] for i in impactor_silicate_particles],
        [i.position_vector[2] for i in impactor_silicate_particles],
        color="red",
        marker="+"
    )
    ax_iron.scatter(
        [i.position_vector[0] for i in impactor_iron_particles],
        [i.position_vector[1] for i in impactor_iron_particles],
        [i.position_vector[2] for i in impactor_iron_particles],
        color="blue",
        marker="+"
    )
    ax_silicate.set_title("Impactor Silicate: Iteration {}".format(time))
    ax_iron.set_title("Impactor Iron: Iteration {}".format(time))
    ax_silicate.set_xlim(-8000e3, 8000e3)
    ax_silicate.set_ylim(-8000e3, 8000e3)
    ax_silicate.set_zlim(-8000e3, 8000e3)
    ax_iron.set_xlim(-8000e3, 8000e3)
    ax_iron.set_ylim(-8000e3, 8000e3)
    ax_iron.set_zlim(-8000e3, 8000e3)
    plt.tight_layout()
    plt.savefig(positions_output_path + "/impactor/{}.png".format(time), format="png")
    plt.close()

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(
        [i.position_vector[0] for i in impactor_silicate_particles] + [i.position_vector[0] for i in
                                                                       impactor_iron_particles],
        [i.position_vector[1] for i in impactor_silicate_particles] + [i.position_vector[1] for i in
                                                                       impactor_iron_particles],
        [i.position_vector[2] for i in impactor_silicate_particles] + [i.position_vector[2] for i in
                                                                       impactor_iron_particles],
        color="red",
        marker="+",
        label="Impactor"
    )
    ax.scatter(
        [i.position_vector[0] for i in target_silicate_particles] + [i.position_vector[0] for i in
                                                                     target_iron_particles],
        [i.position_vector[1] for i in target_silicate_particles] + [i.position_vector[1] for i in
                                                                     target_iron_particles],
        [i.position_vector[2] for i in target_silicate_particles] + [i.position_vector[2] for i in
                                                                     target_iron_particles],
        color="blue",
        marker="+",
        label="Target"
    )
    ax.set_title("Iteration {}".format(time))
    ax.set_xlim(-8000e4, 8000e4)
    ax.set_ylim(-8000e4, 8000e4)
    ax.set_zlim(-8000e4, 8000e4)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(both_output_path + "/{}.png".format(time), format="png")
    plt.close()

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    ax_density = fig.add_subplot(221)
    ax_pressure = fig.add_subplot(222)
    ax_internal_energy = fig.add_subplot(223)
    ax_entropy = fig.add_subplot(224)
    ax_internal_energy.set_xlabel("Distance from Target Center (km)")
    ax_entropy.set_xlabel("Distance from Target Center (km)")
    ax_density.set_title("Density (Iteration {})".format(time))
    ax_pressure.set_title("Pressure (Iteration {})".format(time))
    ax_internal_energy.set_title("Internal Energy (Iteration {})".format(time))
    ax_entropy.set_title("Entropy (Iteration {})".format(time))
    ax_density.grid()
    ax_pressure.grid()
    ax_internal_energy.grid()
    ax_entropy.grid()
    ax_density.scatter(
        target_silicate_distances,
        [i.density for i in target_silicate_particles],
        color="red",
        marker="+",
        label="Silicate"
    )
    ax_density.scatter(
        target_iron_distances,
        [i.density for i in target_iron_particles],
        color="blue",
        marker="+",
        label="Iron"
    )
    # ax_pressure.scatter(
    #     target_silicate_distances,
    #     [i.pressure for i in target_silicate_particles],
    #     color="red",
    #     marker="+",
    #     label="Silicate"
    # )
    # ax_pressure.scatter(
    #     target_iron_distances,
    #     [i.pressure for i in target_iron_particles],
    #     color="blue",
    #     marker="+",
    #     label="Iron"
    # )
    ax_internal_energy.scatter(
        target_silicate_distances,
        [i.internal_energy for i in target_silicate_particles],
        color="red",
        marker="+",
        label="Silicate"
    )
    ax_internal_energy.scatter(
        target_iron_distances,
        [i.internal_energy for i in target_iron_particles],
        color="blue",
        marker="+",
        label="Iron"
    )
    ax_entropy.scatter(
        target_silicate_distances,
        [i.entropy for i in target_silicate_particles],
        color="red",
        marker="+",
        label="Silicate"
    )
    ax_entropy.scatter(
        target_iron_distances,
        [i.entropy for i in target_iron_particles],
        color="blue",
        marker="+",
        label="Iron"
    )
    ax_entropy.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(combined_path + "/target/{}.png".format(time), format="png")
    plt.close()

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    ax_density = fig.add_subplot(221)
    ax_pressure = fig.add_subplot(222)
    ax_internal_energy = fig.add_subplot(223)
    ax_entropy = fig.add_subplot(224)
    ax_internal_energy.set_xlabel("Distance from Target Center (km)")
    ax_entropy.set_xlabel("Distance from Target Center (km)")
    ax_density.set_title("Density (Iteration {})".format(time))
    ax_pressure.set_title("Pressure (Iteration {})".format(time))
    ax_internal_energy.set_title("Internal Energy (Iteration {})".format(time))
    ax_entropy.set_title("Entropy (Iteration {})".format(time))
    ax_density.grid()
    ax_pressure.grid()
    ax_internal_energy.grid()
    ax_entropy.grid()
    ax_density.scatter(
        impactor_silicate_distances,
        [i.density for i in impactor_silicate_particles],
        color="red",
        marker="+",
        label="Silicate"
    )
    ax_density.scatter(
        impactor_iron_distances,
        [i.density for i in impactor_iron_particles],
        color="blue",
        marker="+",
        label="Iron"
    )
    # ax_pressure.scatter(
    #     impactor_silicate_distances,
    #     [i.pressure for i in impactor_silicate_particles],
    #     color="red",
    #     marker="+",
    #     label="Silicate"
    # )
    # ax_pressure.scatter(
    #     impactor_iron_distances,
    #     [i.pressure for i in impactor_iron_particles],
    #     color="blue",
    #     marker="+",
    #     label="Iron"
    # )
    ax_internal_energy.scatter(
        impactor_silicate_distances,
        [i.internal_energy for i in impactor_silicate_particles],
        color="red",
        marker="+",
        label="Silicate"
    )
    ax_internal_energy.scatter(
        impactor_iron_distances,
        [i.internal_energy for i in impactor_iron_particles],
        color="blue",
        marker="+",
        label="Iron"
    )
    ax_entropy.scatter(
        impactor_silicate_distances,
        [i.entropy for i in impactor_silicate_particles],
        color="red",
        marker="+",
        label="Silicate"
    )
    ax_entropy.scatter(
        impactor_iron_distances,
        [i.entropy for i in impactor_iron_particles],
        color="blue",
        marker="+",
        label="Iron"
    )
    ax_entropy.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(combined_path + "/impactor/{}.png".format(time), format="png")
    plt.close()

animate(start_time=start_time, end_time=end_time, interval=interval, output_path=output_path + "/target",
        file_name="target_mode1.mp4")
animate(start_time=start_time, end_time=end_time, interval=interval, output_path=output_path + "/impactor",
        file_name="impactor_mode1.mp4")
animate(start_time=start_time, end_time=end_time, interval=interval, output_path=positions_output_path + "/target",
        file_name="target_profile_mode1.mp4")
animate(start_time=start_time, end_time=end_time, interval=interval, output_path=positions_output_path + "/impactor",
        file_name="impactor_profile_mode1.mp4")
animate(start_time=start_time, end_time=end_time, interval=interval, output_path=both_output_path,
        file_name="both_profile_mode1.mp4")
animate(start_time=start_time, end_time=end_time, interval=interval, output_path=combined_path + "/target",
        file_name="target_combined_profile_mode1.mp4")
animate(start_time=start_time, end_time=end_time, interval=interval, output_path=combined_path + "/impactor",
        file_name="impactor_combined_profile_mode1.mp4")

for i in paths:
    shutil.rmtree(i)
