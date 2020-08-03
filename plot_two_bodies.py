import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from src.identify import ParticleMap

sph_file = "merged_800.dat"
sph_df = pd.read_csv(sph_file, header=None, skiprows=2, delimiter="\t")

pm = ParticleMap(output_path=sph_file, number_expected_bodies=2)
particle_map = pm.solve()

ax = plt.figure().add_subplot(111)
ax.scatter(
    [p.position_vector[0] for p in particle_map if p.assigned_body == "TARGET" and p.label == "PLANET"],
    [p.position_vector[1] for p in particle_map if p.assigned_body == "TARGET" and p.label == "PLANET"],
    marker="+",
    color="blue",
    label="TARGET"
)
ax.scatter(
    [p.position_vector[0] for p in particle_map if p.assigned_body == "IMPACTOR" and p.label == "PLANET"],
    [p.position_vector[1] for p in particle_map if p.assigned_body == "IMPACTOR" and p.label == "PLANET"],
    marker="+",
    color="red",
    label="IMPACTOR"
)
ax.scatter(
    [p.position_vector[0] for p in particle_map if p.label == "DISK"],
    [p.position_vector[1] for p in particle_map if p.label == "DISK"],
    marker="+",
    color="green",
    label="DISK"
)
ax.scatter(
    [p.position_vector[0] for p in particle_map if p.label == "ESCAPE"],
    [p.position_vector[1] for p in particle_map if p.label == "ESCAPE"],
    marker="+",
    color="pink",
    label="ESCAPING"
)
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.grid()
ax.legend()

plt.show()
