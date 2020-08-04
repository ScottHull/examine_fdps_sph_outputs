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
    x=[i.distance / 1000.0 for i in particle_map if i.label == "PLANET"],
    y=[i.entropy for i in particle_map if i.label == "PLANET"],
    tags=[i.particle_id for i in particle_map if i.label == "PLANET"],
    x_label="Distance (km)",
    y_label="Entropy"
)

ax2 = plots.scatter_particles(
    x=[i.position_vector[0] / 1000.0 for i in particle_map],
    y=[i.position_vector[1] / 1000.0 for i in particle_map],
    tags=[i.particle_id for i in particle_map],
    x_label="x (km)",
    y_label="y (km)",
    a=pm.a,
    b=pm.b
)

plt.show()
