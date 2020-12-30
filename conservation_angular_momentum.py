import numpy as np
import matplotlib.pyplot as plt
from src.identify import ParticleMapFromFiles

start_time = 0
end_time = 5000
interval = 200
path = "/scratch/shull4/GI_outfiles"

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("Time Iteration")
ax.set_ylabel("Total Angular Momentum of System")
ax.set_title("System Angular Momentum as Function of Time")
ax.grid()

total_ams = []
for time in np.arange(start_time, end_time + interval, interval):
    particle_map = ParticleMapFromFiles(path=path).read(time=time)
    total_am = sum([p.angular_momentum for p in particle_map])
    total_ams.append(total_am)

ax.plot(
    np.arange(start_time, end_time + interval, interval),
    total_ams,
    linewidth=2.0,
    color="black"
)

plt.savefig("total_angular_momentum.png", format='png')

