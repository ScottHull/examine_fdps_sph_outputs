import pandas as pd
from math import pi
import numpy as np
from src.interpolation import interpolate1d


class Structure:

    def __init__(self, particles, phase="duniteS2"):
        self.particles = particles
        if phase.lower() == "dunites2":
            self.phase_df = pd.read_fwf("src/phase_data/duniteS_vapour_curve.txt", skiprows=1,
                                        names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                               "entropy_sol_liq", "entropy_vap"])

    def calc_vapor_mass_fraction(self, target_label=None):
        num_particles = 0
        vapor_mass_fraction = 0
        for i in self.particles:
            if (target_label is not None and i.label == target_label) or target_label is None:
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
        if num_particles == 0:
            return None
        vapor_mass_fraction = vapor_mass_fraction / num_particles

        return vapor_mass_fraction

    def calc_disk_surface_density(self):
        particles = [x for _, x in sorted(zip([i.distance for i in (self.particles)],
                                              self.particles))]
        surface_densities = []
        for index, p in enumerate(particles):
            if index > 0:
                dM = p.mass - particles[index - 1].mass
                dr = p.distance = particles[index - 1].distance
                sd = dM / (2.0 * pi * p.distance * dr)
                surface_densities.append(sd)
        return surface_densities, [p.distance for p in particles][1:]

    def calc_total_angular_momentum(self, target_label=None):
        angular_momentum = 0
        if (target_label is not None) and target_label == "DISK":
            angular_momentum = sum([np.linalg.norm(p.angular_momentum_vector) for p in self.particles if p.label == "DISK"])
        else:
            angular_momentum = sum([np.linalg.norm(p.angular_momentum_vector) for p in self.particles])
        return angular_momentum

    def calc_total_mass(self, target_label=None):
        mass = 0
        if (target_label is not None) and target_label == "DISK":
            mass = sum(
                [p.mass for p in self.particles if p.label == "DISK"])
        else:
            mass = sum([p.mass for p in self.particles])
        return mass
