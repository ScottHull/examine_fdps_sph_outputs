import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse
from scipy.stats.kde import gaussian_kde
import numpy as np
import moviepy.editor as mpy


def colorcode_orbits(particles, time, a=None, b=None, z=None, center_plot=False):
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)

    escaping_particles = [p for p in particles if p.label == "ESCAPE"]
    disk_particles = [p for p in particles if p.label == "DISK"]
    planet_particles = [p for p in particles if p.label == "PLANET"]

    if z is None:
        ax.scatter(
            [p.position_vector[0] for p in planet_particles],
            [p.position_vector[1] for p in planet_particles],
            c='blue',
            marker="+",
            label="PLANET: {}".format(len(planet_particles))
        )
        ax.scatter(
            [p.position_vector[0] for p in disk_particles],
            [p.position_vector[1] for p in disk_particles],
            c='pink',
            marker="+",
            label="DISK: {}".format(len(disk_particles))
        )
        ax.scatter(
            [p.position_vector[0] for p in escaping_particles],
            [p.position_vector[1] for p in escaping_particles],
            c='red',
            marker="+",
            label="ESCAPE: {}".format(len(escaping_particles))
        )
        if a is not None and b is not None:
            e = Ellipse(xy=(0, 0), width=a * 2.0, height=b * 2.0, alpha=0.3, color="blue")
            ax.add_artist(e)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("PARTICLE ORBITS (iteration: {})".format(time))
        ax.grid()
        ax.legend()

        if center_plot:
            ax.set_xlim(-1e8, 1e8)
            ax.set_ylim(-1e8, 1e8)

    return fig


def eccentricity_plot(particles, time, a):
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)

    escaping_particles = [(p.semi_major_axis / a, p.eccentricity) for p in particles if p.label == "ESCAPE"]
    disk_particles = [(p.semi_major_axis / a, p.eccentricity) for p in particles if p.label == "DISK"]
    planet_particles = [(p.semi_major_axis / a, p.eccentricity) for p in particles if p.label == "PLANET"]

    ax.scatter(
        [p[0] for p in planet_particles],
        [p[1] for p in planet_particles],
        c='blue',
        marker="+",
        label="PLANET: {}".format(len(planet_particles))
    )
    ax.scatter(
        [p[0] for p in disk_particles],
        [p[1] for p in disk_particles],
        c='pink',
        marker="+",
        label="DISK: {}".format(len(disk_particles))
    )
    ax.scatter(
        [p[0] for p in escaping_particles],
        [p[1] for p in escaping_particles],
        c='red',
        marker="+",
        label="ESCAPE: {}".format(len(escaping_particles))
    )
    ax.set_xlabel("a / Rp (semi-major axis / planet radius)")
    ax.set_ylabel("Eccentricity (e)")
    ax.set_title("Iteration: {}".format(time))
    ax.grid()
    ax.legend()

    ax.set_xlim(-1e8, 1e8)
    ax.set_ylim(-1e8, 1e8)

    return fig


def animate(start_time, end_time, interval, path, filename="animation.mp4", fps=15):
    frames = [path + "/{}.png".format(time) for time in np.arange(start_time, end_time + interval, interval)]
    animation = mpy.ImageSequenceClip(frames, fps=fps, load_images=True)
    animation.write_videofile(filename, fps=fps)
