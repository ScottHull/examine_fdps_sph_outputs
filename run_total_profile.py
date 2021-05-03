import os
from src.identify import ParticleMap, ParticleMapFromFiles
from src.combine import CombineFile
from total_profile import controller
import matplotlib.pyplot as plt

time = 2000
number_processes = 100
path_to_outputs = "/scratch/shull4/gi"

combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
f = os.getcwd() + "/merged_{}.dat".format(time)
pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
particle_map = pm.solve()

relative_velocities = controller.get_relative_velocities(particle_map=particle_map)
eccentricities = controller.get_eccentricities(particle_map=particle_map)
periapsis = controller.get_periapsis(particle_map=particle_map)


fig = plt.figure()
ax = fig.add_subplot(111)
ax = controller.plot(ax=ax, ps=relative_velocities)
ax.set_xlabel("Radius")
ax.set_ylabel("Relative Velocity")
ax.set_title("Relative Velocities")
plt.savefig("relative_velocity.png", format='png')

fig = plt.figure()
ax = fig.add_subplot(111)
ax = controller.plot(ax=ax, ps=eccentricities)
ax.set_xlabel("Radius")
ax.set_ylabel("Eccentricity")
ax.set_title("Eccentricities")
plt.savefig("eccentricities.png", format='png')

fig = plt.figure()
ax = fig.add_subplot(111)
ax = controller.plot(ax=ax, ps=periapsis)
ax.set_xlabel("Radius")
ax.set_ylabel("Periapsis")
ax.set_title("Periapsis")
plt.savefig("periapsis.png", format='png')
