import matplotlib.pyplot as plt
import pandas as pd
from math import sqrt

def radius(x_coords, y_coords, z_coords):
    return [sqrt(x**2 + y**2 + z**2) / 1000.0 for x, y, z in zip(x_coords, y_coords, z_coords)]

aneos = [
    (
        "Positive Pressure", "/Users/scotthull/Desktop/merged_25_GI_p_pres.dat"
    ),
    (
        "Positive Pressure & Positive Energy", "/Users/scotthull/Desktop/merged_25_GI_p_eng_p_pres.dat"
    )
]

tillotson = [
    (
        "No Positive Pressure/Energy Rule", "/Users/scotthull/Desktop/merged_25_tillotson2.dat"
    ),
    (
        "Positive Pressure/Energy Rule", "/Users/scotthull/Desktop/merged_25_tillotson.dat"
    )
]

fig = plt.figure()
fig2 = plt.figure()
ax_aneos = fig.add_subplot(111)
ax_tillotson = fig2.add_subplot(111)

plotting_index = 10
plotting_label = "Internal Energy"

for i in aneos:
    df = pd.read_csv(i[1], sep='\t', skiprows=2, header=None)
    r = radius(x_coords=df[3], y_coords=df[4], z_coords=df[5])
    ax_aneos.scatter(r, df[plotting_index], marker="+", label=i[0])
for i in tillotson:
    df = pd.read_csv(i[1], sep='\t', skiprows=2, header=None)
    r = radius(x_coords=df[3], y_coords=df[4], z_coords=df[5])
    ax_tillotson.scatter(r, df[plotting_index], marker="+", label=i[0])

ax_aneos.set_xlabel("Radius (km)")
ax_tillotson.set_xlabel("Radius (km)")
ax_aneos.set_ylabel(plotting_label)
ax_tillotson.set_ylabel(plotting_label)
ax_aneos.set_title("ANEOS")
ax_tillotson.set_title("Tillotson")
ax_aneos.grid()
ax_tillotson.grid()
ax_aneos.legend()
ax_tillotson.legend()

plt.tight_layout()
plt.show()
