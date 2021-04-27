from src.orbitalelements import Particle
from src.orbitalelements_3d import Particle as Particle3D


position_vector = [-1720556620, -1113064066, -181929491.6]
velocity_vector = [-8534.470778, -5569.50806, -886.787378]

p = Particle(
    particle_name=0,
    particle_id=0,
    position_vector=position_vector,
    velocity_vector=velocity_vector,
    mass=5.59E+19,
    density=5,
    internal_energy=18592113.85,
    mass_grav_body=5e23,
    entropy=5793.538662,
    temperature=2000
)
p.pressure = 10

p3d = Particle3D(
    particle_name=0,
    particle_id=0,
    position_vector=position_vector,
    velocity_vector=velocity_vector,
    mass=5.59E+19,
    density=5,
    internal_energy=18592113.85,
    mass_grav_body=5e23,
    entropy=5793.538662,
    temperature=2000
)
p.pressure = 10

