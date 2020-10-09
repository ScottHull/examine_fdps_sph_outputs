import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats.kde import gaussian_kde
import numpy as np


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


def scatter_particles(x, y, tags, x_label, y_label, a=None, b=None):
    ax = plt.figure().add_subplot(111)

    target_silicate_x = [x[index] for index, i in enumerate(tags) if i == 0]
    target_silicate_y = [y[index] for index, i in enumerate(tags) if i == 0]
    target_iron_x = [x[index] for index, i in enumerate(tags) if i == 1]
    target_iron_y = [y[index] for index, i in enumerate(tags) if i == 1]
    impactor_silicate_x = [x[index] for index, i in enumerate(tags) if i == 2]
    impactor_silicate_y = [y[index] for index, i in enumerate(tags) if i == 2]
    impactor_iron_x = [x[index] for index, i in enumerate(tags) if i == 3]
    impactor_iron_y = [y[index] for index, i in enumerate(tags) if i == 3]

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

    return ax


def colorcode_orbits(particles, a, b):
    ax = plt.figure().add_subplot(111)
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
        [p.position_vector[0] for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        [p.position_vector[1] for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
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

    return ax


def plot_eccentricities(particles, a, b):
    ax = plt.figure().add_subplot(111)
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
    ax.legend()

    return ax


def plot_eccentricity_elements(particles, a, b):
    fig = plt.figure()
    ax_angular_momentum = fig.add_subplot(411)
    ax_orbital_energy = fig.add_subplot(412)
    ax_alpha = fig.add_subplot(413)
    ax_mass_reduced = fig.add_subplot(414)
    ax_angular_momentum.set_ylabel("L")
    ax_orbital_energy.set_ylabel("E")
    ax_alpha.set_ylabel("alpha")
    ax_mass_reduced.set_ylabel("mass_red")
    ax_angular_momentum.grid()
    ax_orbital_energy.grid()
    ax_alpha.grid()
    ax_mass_reduced.grid()

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

    ax_alpha.scatter(
        [p.distance for p in particles if p.eccentricity > 1.0],
        [p.alpha for p in particles if p.eccentricity > 1.0],
        c='red',
        marker="+",
        label="ESCAPE"
    )
    ax_alpha.scatter(
        [p.distance for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        [p.alpha for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        c='blue',
        marker="+",
        label="DISTANCE WITHIN PLANET"
    )
    ax_alpha.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        [p.alpha for p in particles if
         p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        c='green',
        marker="+",
        label="PERIAPSIS INSIDE PLANET"
    )
    ax_alpha.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        [p.alpha for p in particles if
         p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        c='pink',
        marker="+",
        label="DISK"
    )

    ax_mass_reduced.scatter(
        [p.distance for p in particles if p.eccentricity > 1.0],
        [p.mass_reduced for p in particles if p.eccentricity > 1.0],
        c='red',
        marker="+",
        label="ESCAPE"
    )
    ax_mass_reduced.scatter(
        [p.distance for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        [p.mass_reduced for p in particles if
         abs(p.position_vector[0]) <= a and abs(p.position_vector[2]) <= a and abs(p.position_vector[1]) <= b],
        c='blue',
        marker="+",
        label="DISTANCE WITHIN PLANET"
    )
    ax_mass_reduced.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        [p.mass_reduced for p in particles if
         p.eccentricity <= 1.0 and abs(p.periapsis) <= a and p.distance > a],
        c='green',
        marker="+",
        label="PERIAPSIS INSIDE PLANET"
    )
    ax_mass_reduced.scatter(
        [p.distance for p in particles if p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        [p.mass_reduced for p in particles if
         p.eccentricity <= 1.0 and abs(p.periapsis) > a],
        c='pink',
        marker="+",
        label="DISK"
    )

    ax_angular_momentum.legend()
    ax_mass_reduced.set_xlabel("Distance")
    return fig


def plot_velocity(particles, a, b):
    ax = plt.figure().add_subplot(111)
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

    return ax

def plot_energies(particles, a, b):
    ax = plt.figure().add_subplot(111)
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

    return ax
