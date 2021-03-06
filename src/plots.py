import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse
from scipy.stats.kde import gaussian_kde
import numpy as np
import moviepy.editor as mpy


def plot_heatmap(x, y, z, a, b, center):
    ax = plt.figure().add_subplot(111)
    sc = ax.scatter(x, y, c=z, marker="+")
    e = Ellipse(xy=(0, 0), width=a * 2.0, height=b * 2.0, alpha=0.2, color="blue")
    ax.scatter(center[0], center[1], marker="*", color="red", s=200)
    ax.add_artist(e)
    cbar = plt.colorbar(sc, ax=ax)

    return ax


def plot_particle_density_heatmap(x, y):
    x = np.array(x)
    y = np.array(y)
    k = gaussian_kde(np.vstack([x, y]))
    xi, yi = np.mgrid[x.min():x.max():x.size ** 0.5 * 1j, y.min():y.max():y.size ** 0.5 * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    ax = plt.figure().add_subplot(111)
    # alpha=0.5 will make the plots semitransparent
    ax.pcolormesh(xi, yi, zi.reshape(xi.shape), alpha=0.5)

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())

    return ax


def scatter_particles(x, y, tags, x_label, y_label, z=None, z_label=None, a=None, b=None, center_plot=False):
    target_silicate_x = [x[index] for index, i in enumerate(tags) if i == 0]
    target_silicate_y = [y[index] for index, i in enumerate(tags) if i == 0]
    target_silicate_z = [z[index] for index, i in enumerate(tags) if i == 0]
    target_iron_x = [x[index] for index, i in enumerate(tags) if i == 1]
    target_iron_y = [y[index] for index, i in enumerate(tags) if i == 1]
    target_iron_z = [z[index] for index, i in enumerate(tags) if i == 1]
    impactor_silicate_x = [x[index] for index, i in enumerate(tags) if i == 2]
    impactor_silicate_y = [y[index] for index, i in enumerate(tags) if i == 2]
    impactor_silicate_z = [z[index] for index, i in enumerate(tags) if i == 2]
    impactor_iron_x = [x[index] for index, i in enumerate(tags) if i == 3]
    impactor_iron_y = [y[index] for index, i in enumerate(tags) if i == 3]
    impactor_iron_z = [z[index] for index, i in enumerate(tags) if i == 3]

    fig = plt.figure()

    if z is None:
        ax = fig.add_subplot(111)
        ax.scatter(target_silicate_x, target_silicate_y, marker="+", color="red", label="Target Silicate")
        ax.scatter(target_iron_x, target_iron_y, marker="+", color="blue", label="Target Iron")
        ax.scatter(impactor_silicate_x, impactor_silicate_y, marker="+", color="green", label="Impactor Silicate")
        ax.scatter(impactor_iron_x, impactor_iron_y, marker="+", color="pink", label="Impactor Iron")

        if a is not None and b is not None:
            e = Ellipse(xy=(0, 0), width=a * 2.0, height=b * 2.0, alpha=0.3, color="blue")
            ax.add_artist(e)

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.grid()
        ax.legend()

        if center_plot:
            ax.set_xlim(-2e8, 2e8)
            ax.set_ylim(-2e8, 2e8)

    else:
        ax = Axes3D(fig)
        ax.scatter(target_silicate_x, target_silicate_y, target_silicate_z, marker="+", color="red",
                   label="Target Silicate")
        ax.scatter(target_iron_x, target_iron_y, target_iron_z, marker="+", color="blue", label="Target Iron")
        ax.scatter(impactor_silicate_x, impactor_silicate_y, impactor_silicate_z, marker="+", color="green",
                   label="Impactor Silicate")
        ax.scatter(impactor_iron_x, impactor_iron_y, impactor_iron_z, marker="+", color="pink", label="Impactor Iron")

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_zlabel(z_label)
        ax.grid()
        ax.legend(loc="lower right")

        if center_plot:
            ax.set_xlim(-2e8, 2e8)
            ax.set_ylim(-2e8, 2e8)
            ax.set_zlim(-2e8, 2e8)

    return fig


def colorcode_orbits(particles, a, b, z=None, center_plot=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if z is None:
        ax.scatter(
            [p.position_vector[0] for p in particles if p.eccentricity > 1.0],
            [p.position_vector[1] for p in particles if p.eccentricity > 1.0],
            c='red',
            marker="+",
            label="ESCAPE"
        )
        ax.scatter(
            [p.position_vector[0] for p in particles if
             abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
            [p.position_vector[1] for p in particles if
             abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
            c='blue',
            marker="+",
            label="DISTANCE WITHIN PLANET"
        )
        ax.scatter(
            [p.position_vector[0] for p in particles if
             p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
            [p.position_vector[1] for p in particles if
             p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
            c='green',
            marker="+",
            label="PERIAPSIS INSIDE PLANET"
        )
        ax.scatter(
            [p.position_vector[0] for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
            [p.position_vector[1] for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
            c='pink',
            marker="+",
            label="DISK"
        )
        if a is not None and b is not None:
            e = Ellipse(xy=(0, 0), width=a * 2.0, height=b * 2.0, alpha=0.3, color="blue")
            ax.add_artist(e)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("PARTICLE ORBITS")
        ax.grid()
        ax.legend()

        if center_plot:
            ax.set_xlim(-2e8, 2e8)
            ax.set_ylim(-2e8, 2e8)

    else:
        ax = Axes3D(fig)
        ax.scatter(
            [p.position_vector[0] for p in particles if p.eccentricity > 1.0],
            [p.position_vector[1] for p in particles if p.eccentricity > 1.0],
            [p.position_vector[2] for p in particles if p.eccentricity > 1.0],
            c='red',
            marker="+",
            label="ESCAPE"
        )
        ax.scatter(
            [p.position_vector[0] for p in particles if
             abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
            [p.position_vector[1] for p in particles if
             abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
            [p.position_vector[2] for p in particles if
             abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
            c='blue',
            marker="+",
            label="DISTANCE WITHIN PLANET"
        )
        ax.scatter(
            [p.position_vector[0] for p in particles if
             p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
            [p.position_vector[1] for p in particles if
             p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
            [p.position_vector[2] for p in particles if
             p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
            c='green',
            marker="+",
            label="PERIAPSIS INSIDE PLANET"
        )
        ax.scatter(
            [p.position_vector[0] for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
            [p.position_vector[1] for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
            [p.position_vector[2] for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
            c='pink',
            marker="+",
            label="DISK"
        )

        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_title("PARTICLE ORBITS")
        ax.grid()
        ax.legend(loc="lower right")

        if center_plot:
            ax.set_xlim(-2e8, 2e8)
            ax.set_ylim(-2e8, 2e8)
            ax.set_zlim(-2e8, 2e8)

    return fig


def plot_eccentricities(particles, a, b):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(
        [p.distance for p in particles if p.eccentricity > 1.0],
        [p.eccentricity for p in particles if p.eccentricity > 1.0],
        c='red',
        marker="+",
        label="ESCAPE"
    )
    ax.scatter(
        [p.distance for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        [p.eccentricity for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        c='blue',
        marker="+",
        label="DISTANCE WITHIN PLANET"
    )
    ax.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        [p.eccentricity for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        c='green',
        marker="+",
        label="PERIAPSIS INSIDE PLANET"
    )
    ax.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        [p.eccentricity for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        c='pink',
        marker="+",
        label="DISK"
    )

    ax.set_xlabel("Radial Distance From Earth Center")
    ax.set_ylabel("Orbital Eccentricity")
    ax.set_title("PARTICLE ORBITAL ECCENTRICITIES")
    ax.grid()
    ax.legend(loc="upper right")

    return fig


def plot_eccentricity_elements(particles, a, b):
    fig = plt.figure()
    ax_angular_momentum = fig.add_subplot(211)
    ax_orbital_energy = fig.add_subplot(212)
    ax_angular_momentum.set_ylabel("L")
    ax_orbital_energy.set_ylabel("E")
    ax_angular_momentum.grid()
    ax_orbital_energy.grid()

    ax_angular_momentum.scatter(
        [p.distance for p in particles if p.eccentricity > 1.0],
        [np.linalg.norm(p.angular_momentum_vector) for p in particles if p.eccentricity > 1.0],
        c='red',
        marker="+",
        label="ESCAPE"
    )
    ax_angular_momentum.scatter(
        [p.distance for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        [np.linalg.norm(p.angular_momentum_vector) for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        c='blue',
        marker="+",
        label="DISTANCE WITHIN PLANET"
    )
    ax_angular_momentum.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        [np.linalg.norm(p.angular_momentum_vector) for p in particles if
         p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        c='green',
        marker="+",
        label="PERIAPSIS INSIDE PLANET"
    )
    ax_angular_momentum.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        [np.linalg.norm(p.angular_momentum_vector) for p in particles if
         p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        c='pink',
        marker="+",
        label="DISK"
    )

    ax_orbital_energy.scatter(
        [p.distance for p in particles if p.eccentricity > 1.0],
        [p.orbital_energy for p in particles if p.eccentricity > 1.0],
        c='red',
        marker="+",
        label="ESCAPE"
    )
    ax_orbital_energy.scatter(
        [p.distance for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        [p.orbital_energy for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        c='blue',
        marker="+",
        label="DISTANCE WITHIN PLANET"
    )
    ax_orbital_energy.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        [p.orbital_energy for p in particles if
         p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        c='green',
        marker="+",
        label="PERIAPSIS INSIDE PLANET"
    )
    ax_orbital_energy.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        [p.orbital_energy for p in particles if
         p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        c='pink',
        marker="+",
        label="DISK"
    )

    ax_angular_momentum.legend()
    ax_orbital_energy.set_xlabel("Distance")
    return fig


def plot_velocity(particles, a, b):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(
        [p.distance for p in particles if p.eccentricity > 1.0],
        [np.linalg.norm(p.velocity_vector) for p in particles if p.eccentricity > 1.0],
        c='red',
        marker="+",
        label="ESCAPE"
    )
    ax.scatter(
        [p.distance for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        [np.linalg.norm(p.velocity_vector) for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        c='blue',
        marker="+",
        label="DISTANCE WITHIN PLANET"
    )
    ax.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        [np.linalg.norm(p.velocity_vector) for p in particles if
         p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        c='green',
        marker="+",
        label="PERIAPSIS INSIDE PLANET"
    )
    ax.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        [np.linalg.norm(p.velocity_vector) for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        c='pink',
        marker="+",
        label="DISK"
    )

    ax.set_xlabel("Radial Distance From Earth Center")
    ax.set_ylabel("Orbital Velocity")
    ax.set_title("PARTICLE ORBITAL VELOCITIES")
    ax.grid()
    ax.legend()

    return fig


def plot_energies(particles, a, b):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(
        [p.distance for p in particles if p.eccentricity > 1.0],
        [p.potential_energy for p in particles if p.eccentricity > 1.0],
        c='red',
        marker="+",
        label="ESCAPE (PE)"
    )
    ax.scatter(
        [p.distance for p in particles if p.eccentricity > 1.0],
        [p.kinetic_energy for p in particles if p.eccentricity > 1.0],
        c='red',
        marker="o",
        label="ESCAPE (KE)"
    )
    ax.scatter(
        [p.distance for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        [p.potential_energy for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        c='blue',
        marker="+",
        label="DISTANCE WITHIN PLANET (PE)"
    )
    ax.scatter(
        [p.distance for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        [p.kinetic_energy for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        c='blue',
        marker="o",
        label="DISTANCE WITHIN PLANET (KE)"
    )
    ax.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        [p.potential_energy for p in particles if
         p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        c='green',
        marker="+",
        label="PERIAPSIS INSIDE PLANET (PE)"
    )
    ax.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        [p.kinetic_energy for p in particles if
         p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        c='green',
        marker="o",
        label="PERIAPSIS INSIDE PLANET (KE)"
    )
    ax.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        [p.potential_energy for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        c='pink',
        marker="+",
        label="DISK (PE)"
    )
    ax.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        [p.kinetic_energy for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        c='pink',
        marker="o",
        label="DISK (KE)"
    )

    ax.set_xlabel("Radial Distance From Earth Center")
    ax.set_ylabel("Orbital Energies")
    ax.set_title("PARTICLE ORBITAL ENERGIES")
    ax.grid()
    ax.legend()

    return fig


def animate(start_time, end_time, interval, path, filename="animation.mp4", fps=30):
    frames = [path + "/{}.png".format(time) for time in np.arange(start_time, end_time + interval, interval)]
    animation = mpy.ImageSequenceClip(frames, fps=fps, load_images=True)
    animation.write_videofile(filename, fps=fps)


def plot_vfm(phase_curve_1_x, phase_curve_1_y, particles_color, phase_curve_2_x, phase_curve_2_y, particles_x,
             particles_y, xlabel, ylabel, cbar_label, phase_curve_1_label="sol-liq", phase_curve_2_label="liq-vap"):
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    cm = plt.cm.get_cmap('RdYlBu')

    sc = ax.scatter(particles_x, particles_y, c=particles_color, marker="+", cmap=cm)
    cbar = plt.colorbar(sc)
    ax.plot(phase_curve_1_x, phase_curve_1_y, linewidth=2.0, label=phase_curve_1_label)
    ax.plot(phase_curve_2_x, phase_curve_2_y, linewidth=2.0, label=phase_curve_2_label)
    cbar.set_label(cbar_label)
    ax.set_title("VAPOR MASS FRACTION")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    ax.legend()
    # fig.tight_layout()

    return fig
