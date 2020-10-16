import matplotlib.pyplot as plt
import pandas as pd


paths = [
    ("src/eos/dunite.rho_u.txt", "dunite"),
    ("src/eos/dunite3.rho_u.txt", "dunite3"),
    ("src/eos/duniteS2.rho_u.txt", "duniteS2")
]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("Density")
ax.set_ylabel("Internal Energy")
ax.set_title("Compare EOS")

for i in paths:
    path, name = i
    eos_df = pd.read_fwf(path, header=None, skiprows=2)
    density, internal_energy = eos_df[0], eos_df[1]
    ax.scatter(
        density,
        internal_energy,
        marker="+",
        label=name
    )

ax.grid()
ax.legend()

plt.show()
