import os
from src.interpolation import interpolate1d, NearestNeighbor1D
from src.identify import ParticleMap
from src.combine import CombineFile
from src.structure import Structure
import matplotlib.pyplot as plt
import shutil
from random import randint

time = 5000
number_processes = 100
path_to_outputs = "/scratch/shull4/GI"
max_plot_particles = 10

combined_file = CombineFile(num_processes=number_processes, time=time, output_path=path_to_outputs).combine()
f = os.getcwd() + "/merged_{}.dat".format(time)
pm = ParticleMap(output_path=f, center_on_target_iron=True, plot=False, relative_velocity=True, center_plot=True)
particle_map = pm.solve()
s = Structure(particles=particle_map, phase="duniteS2")


def get_phase_curve_temperature_from_entropy(entropy_value, entropy_list, temperature_list):
    index = list(entropy_list).index(entropy_value)
    return list(temperature_list)[index]


rand_particles = []
sol_liq_interp = []
liq_vap_interp = []

selected_particles = 0
while selected_particles < max_plot_particles:
    rand = randint(0, len(particle_map))
    p = particle_map[rand]
    if p.label == "DISK" and p.particle_id % 2 == 0:
        rand_particles.append(p)
        # entropy_liq = interpolate1d(val=p.temperature, val_array=s.phase_df['temperature'],
        #                             interp_array=s.phase_df['entropy_sol_liq'])
        # entropy_vap = interpolate1d(val=p.temperature, val_array=s.phase_df['temperature'],
        #                             interp_array=s.phase_df['entropy_vap'])
        # temp_liq = get_phase_curve_temperature_from_entropy(entropy_value=entropy_liq,
        #                                                     temperature_list=s.phase_df['temperature'],
        #                                                     entropy_list=s.phase_df['entropy_sol_liq'])
        # temp_vap = get_phase_curve_temperature_from_entropy(entropy_value=entropy_liq,
        #                                                     temperature_list=s.phase_df['temperature'],
        #                                                     entropy_list=s.phase_df['entropy_vap'])
        nearest_temp_index = NearestNeighbor1D().neighbor_index(given_val=p.temperature,
                                                                array=list(s.phase_df['temperature']))
        entropy_liq = s.phase_df['entropy_sol_liq'][nearest_temp_index]
        entropy_vap = s.phase_df['entropy_vap'][nearest_temp_index]
        sol_liq_interp.append((entropy_liq, s.phase_df['temperature'][nearest_temp_index]))
        liq_vap_interp.append((entropy_vap, s.phase_df['temperature'][nearest_temp_index]))
        selected_particles += 1

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
cm = plt.cm.get_cmap('RdYlBu')
sc = ax.scatter(
    [p.entropy for p in particle_map if p.label == "DISK" and p.particle_id % 2 == 0],
    [p.temperature for p in particle_map if p.label == "DISK" and p.particle_id % 2 == 0],
    c=[p.distance for p in particle_map if p.label == "DISK" and p.particle_id % 2 == 0],
    marker="+",
    alpha=0.1
)
cbar = plt.colorbar(sc)
ax.scatter(
    [p.entropy for p in rand_particles],
    [p.temperature for p in rand_particles],
    color="black",
    marker="*",
    alpha=1,
)
for index, p in enumerate(rand_particles):
    liq = sol_liq_interp[index]
    vap = liq_vap_interp[index]
    ax.plot(
        [p.entropy, liq[0]],
        [p.temperature, liq[1]],
        color='red',
        linewidth=1.5,
    )
    ax.plot(
        [p.entropy, vap[0]],
        [p.temperature, vap[1]],
        color='purple',
        linewidth=1.5
    )

ax.plot(
    s.phase_df['entropy_sol_liq'],
    s.phase_df['temperature'],
    linewidth=2.0,
    label="Sol-Liq Phase Boundary"
)
ax.plot(
    s.phase_df['entropy_vap'],
    s.phase_df['temperature'],
    linewidth=2.0,
    label="Liq-Vap Phase Boundary"
)
ax.set_xlim(0, 10000)
cbar.set_label("Radial Distance from Target Center")
ax.set_xlabel("Entropy")
ax.set_ylabel("Temperature")
ax.grid()
ax.legend()
fig.savefig("verify_entropy_interp.png", format="png")
