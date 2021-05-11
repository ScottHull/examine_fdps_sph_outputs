import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("/Users/scotthull/Documents/examine_fdps_sph_outputs/results.00000_00001_00000.dat",
                 sep='\t', header=None)
id = df[1]
x = df[3]
y = df[4]
z = df[5]

mantle_x = []
mantle_y = []
mantle_z = []
core_x = []
core_y = []
core_z = []

for index, i in enumerate(id):
    if i == 0:
        mantle_x.append(x[index])
        mantle_y.append(y[index])
        mantle_z.append(z[index])
    elif i == 1:
        core_x.append(x[index])
        core_y.append(y[index])
        core_z.append(z[index])

fig = plt.figure()
ax = fig.add_subplot(111)
for index, i in enumerate(mantle_x):
    if -3000000 < i < 3000000:
        ax.scatter(y[index], z[index], color='black')

plt.show()
