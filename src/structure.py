import pandas as pd
from math import pi
import numpy as np
from src.interpolation import interpolate1d


class Structure:

    def __init__(self, disk_particles, phase="duniteS2"):
        self.disk_particles = disk_particles
        if phase.lower() == "dunites2":
            self.phase_df = pd.read_fwf("src/phase_data/duniteS_vapour_curve.txt", skiprows=1,
                                        names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                               "entropy_sol_liq", "entropy_vap"])

    def calc_vapor_mass_fraction(self):
        num_particles = 0
        vapor_mass_fraction = 0
        for i in self.disk_particles:
            num_particles += 1
            entropy_i = i.entropy
            temperature_i = i.temperature
            entropy_liq = interpolate1d(val=temperature_i, val_array=self.phase_df['temperature'],
                                        interp_array=self.phase_df['entropy_sol_liq'])
            entropy_vap = interpolate1d(val=temperature_i, val_array=self.phase_df['temperature'],
                                        interp_array=self.phase_df['entropy_vap'])
            if entropy_i < entropy_liq:
                vapor_mass_fraction += 0.0
            elif entropy_liq <= entropy_i <= entropy_vap:
                vapor_mass_fraction += (entropy_i - entropy_liq) / (entropy_vap - entropy_liq)
            elif entropy_i > entropy_vap:
                vapor_mass_fraction += 1.0
        vapor_mass_fraction = vapor_mass_fraction / num_particles

        return vapor_mass_fraction

    def calc_disk_surface_density(self):
        particles = [x for _, x in sorted(zip([i.distance for i in self.disk_particles],
                                                                    self.disk_particles))]
        surface_densities = []
        for index, p in enumerate(particles):
            if index > 0:
                dM = p.mass - particles[index - 1].mass
                dr = p.distance = particles[index - 1].distance
                sd = dM / (2.0 * pi * p.distance * dr)
                surface_densities.append(sd)
        return surface_densities, [p.distance for p in particles][1:]



