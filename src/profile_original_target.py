import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from src.identify import ParticleMap
import src.plots as plots

sph_file = "merged_1000.dat"
sph_df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")

pm = ParticleMap(output_path=sph_file, number_expected_bodies=1)
particle_map = pm.solve()

ax1 = plots.scatter_particles(
    x=[i.distance / 1000.0 for i in particle_map if i.particle_id == 0 or i.particle_id == 1],
    y=[i.density for i in particle_map if i.particle_id == 0 or i.particle_id == 1],
    tags=[i.particle_id for i in particle_map if i.particle_id == 0 or i.particle_id == 1],
    x_label="Distance (km)",
    y_label="Density"
)

ax2 = plots.scatter_particles(
    x=[i.distance / 1000.0 for i in particle_map if i.particle_id == 0 or i.particle_id == 1],
    y=[i.internal_energy for i in particle_map if i.particle_id == 0 or i.particle_id == 1],
    tags=[i.particle_id for i in particle_map if i.particle_id == 0 or i.particle_id == 1],
    x_label="Distance (km)",
    y_label="Internal Energy"
)

ax3 = plots.scatter_particles(
    x=[i.position_vector[0] / 1000.0 for i in particle_map],
    y=[i.position_vector[1] / 1000.0 for i in particle_map],
    tags=[i.particle_id for i in particle_map],
    x_label="x (km)",
    y_label="y (km)",
    a=pm.a / 1000.0,
    b=pm.b / 1000.0
)

plt.show()
