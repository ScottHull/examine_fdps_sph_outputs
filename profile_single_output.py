import matplotlib.pyplot as plt
from src.identify import ParticleMapFromFiles

path = "/Users/scotthull/Desktop"
time = 5000

particle_map = ParticleMapFromFiles(path=path).read(time=time)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(
    [i.distance / 1000.0 for i in particle_map if i.label == "ESCAPE"],
    [i.entropy for i in particle_map if i.label == "ESCAPE"],
    color="red",
    marker="+",
    label="Escaping",
)
ax.scatter(
    [i.distance / 1000.0 for i in particle_map if i.label == "DISK"],
    [i.entropy for i in particle_map if i.label == "DISK"],
    color="pink",
    marker="+",
    label="Disk",
)
ax.scatter(
    [i.distance / 1000.0 for i in particle_map if i.label == "PLANET"],
    [i.entropy for i in particle_map if i.label == "PLANET"],
    color="blue",
    marker="+",
    label="Planet",
)
ax.set_xlabel("Distance from Target Center (km)")
ax.set_ylabel("Entropy")
ax.set_title("Distance vs. Entropy at Iteration: {}".format(time))
ax.grid()
ax.legend()
plt.show()
