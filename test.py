import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

path = "/Users/scotthull/Desktop/2000.csv"
df = pd.read_csv(path)

planet = []
disk = []
escape = []
ids = []

for row in df.index:
    label = df['label'][row]
    position = [df['x'][row], df['y'][row], df['z'][row]]
    i = df['particle_id'][row]
    ids.append(i)
    if label == "PLANET":
        planet.append(position)
    elif label == "DISK":
        disk.append(position)
    else:
        escape.append(position)

print([item for item, count in Counter(ids).items() if count > 1])
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    [p[0] for p in planet],
    [p[1] for p in planet],
    marker="+",
    color="blue",
    label="PLANET"
)
ax.scatter(
    [p[0] for p in disk],
    [p[1] for p in disk],
    marker="+",
    color="pink",
    label="DISK"
)
ax.scatter(
    [p[0] for p in escape],
    [p[1] for p in escape],
    marker="+",
    color="red",
    label="ESCAPE"
)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim(-1e8, 1e8)
ax.set_ylim(-1e8, 1e8)
ax.grid()
ax.legend(loc='upper right')
plt.show()
