import numpy as np


def particle_location(particle_map):
    planet = [p for p in particle_map if p.label == "PLANET"]
    disk = [p for p in particle_map if p.label == "DISK"]
    escape = [p for p in particle_map if p.label == "ESCAPE"]

    return {
        'planet': planet,
        'disk': disk,
        'escape': escape
    }


def __relative_velocity(p):
    return np.linalg.norm([p.relative_velocity_vector[0], p.relative_velocity_vector[1], p.relative_velocity_vector[2]])


def get_relative_velocities(particle_map):
    pm = particle_location(particle_map=particle_map)
    return {
        'planet': [(p.distance, __relative_velocity(p)) for p in pm['planet']],
        'disk': [(p.distance, __relative_velocity(p)) for p in pm['disk']],
        'escape': [(p.distance, __relative_velocity(p)) for p in pm['escape']],
    }


def get_eccentricities(particle_map):
    pm = particle_location(particle_map=particle_map)
    return {
        'planet': [(p.distance, p.eccentricity) for p in pm['planet']],
        'disk': [(p.distance, p.eccentricity) for p in pm['disk']],
        'escape': [(p.distance, p.eccentricity) for p in pm['escape']],
    }

def get_periapsis(particle_map):
    pm = particle_location(particle_map=particle_map)
    return {
        'planet': [(p.distance, p.periapsis) for p in pm['planet']],
        'disk': [(p.distance, p.periapsis) for p in pm['disk']],
        'escape': [(p.distance, p.periapsis) for p in pm['escape']],
    }

def plot(ax, ps):
    ax.scatter(
        [p[0] for p in ps['planet']],
        [p[1] for p in ps['planet']],
        marker="+",
        color='blue',
        label='PLANET ({})'.format(len(ps['planet']))
    )
    ax.scatter(
        [p[0] for p in ps['disk']],
        [p[1] for p in ps['disk']],
        marker="+",
        color='pink',
        label='DISK ({})'.format(len(ps['disk']))
    )
    ax.scatter(
        [p[0] for p in ps['escape']],
        [p[1] for p in ps['escape']],
        marker="+",
        color='red',
        label='ESCAPE ({})'.format(len(ps['disk']))
    )
    ax.grid()
    ax.legend()
    return ax
