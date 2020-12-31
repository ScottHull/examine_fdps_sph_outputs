import numpy as np
import matplotlib.pyplot as plt
from src.identify import ParticleMapFromFiles

start_time = 0
end_time = 5000
interval = 200
path = "/scratch/shull4/GI_outfiles"

fig = plt.figure()
fig.set_size_inches(18.5, 10.5)
ax = fig.add_subplot(111)
ax.set_xlabel("Time Iteration")
ax.set_ylabel("Total Angular Momentum of System")
ax.set_title("System Angular Momentum as Function of Time")
ax.grid()

total_ams = []
total_momentum_x = []
total_momentum_y = []
total_momentum_z = []
total_momentums = []
for time in np.arange(start_time, end_time + interval, interval):
    particle_map = ParticleMapFromFiles(path=path).read(time=time)
    total_am = sum([p.angular_momentum for p in particle_map])
    momentum_x = sum([p.momentum_vector[0] for p in particle_map])
    momentum_y = sum([p.momentum_vector[1] for p in particle_map])
    momentum_z = sum([p.momentum_vector[2] for p in particle_map])
    total_momentum = momentum_x + momentum_y + momentum_z
    total_mass = sum([p.mass for p in particle_map])
    total_ams.append(total_am / total_mass)
    total_momentum_x.append(momentum_x / total_mass)
    total_momentum_y.append(momentum_y / total_mass)
    total_momentum_z.append(momentum_z / total_mass)
    total_momentums.append(total_momentum / total_mass)

ax.plot(
    np.arange(start_time, end_time + interval, interval),
    total_ams,
    linewidth=2.0,
    color="black"
)

plt.savefig("total_angular_momentum.png", format='png')
fig.clear()

fig = plt.figure()
fig.set_size_inches(18.5, 10.5)
ax_x = fig.add_subplot(311)
ax_y = fig.add_subplot(312)
ax_z = fig.add_subplot(313)
ax_z.set_xlabel("Time Iteration")
ax_x.set_ylabel("Normalized Momentum")
ax_y.set_ylabel("Normalized Momentum")
ax_z.set_ylabel("Normalized Momentum")
ax_x.set_title("x-momentum")
ax_y.set_title("y-momentum")
ax_z.set_title("z-momentum")
ax_x.grid()
ax_y.grid()
ax_z.grid()
ax_x.plot(
    np.arange(start_time, end_time + interval, interval),
    total_momentum_x,
    linewidth=2.0,
    label="x"
)
ax_y.plot(
    np.arange(start_time, end_time + interval, interval),
    total_momentum_y,
    linewidth=2.0,
    label="y"
)
ax_z.plot(
    np.arange(start_time, end_time + interval, interval),
    total_momentum_z,
    linewidth=2.0,
    label="z"
)
plt.savefig("total_momentum_components.png", format='png')
fig.clear()

fig = plt.figure()
fig.set_size_inches(18.5, 10.5)
ax = fig.add_subplot(111)
ax.set_xlabel("Time Iteration")
ax.set_ylabel("Total Normalized Momentum of System")
ax.set_title("System Total Momentum as Function of Time")
ax.grid()
ax.plot(
    np.arange(start_time, end_time + interval, interval),
    total_momentums,
    linewidth=2.0,
)
plt.savefig("total_momentum.png", format='png')
fig.clear()

