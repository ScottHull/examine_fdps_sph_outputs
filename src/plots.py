import matplotlib.pyplot as plt
import numpy as np

def plot_heatmap(x, y, z):
    ax = plt.figure().add_subplot(111)
    # y, x = np.meshgrid(x, y)
    z = z[:-1]
    c = ax.pcolormesh(np.array(x), np.array(y), np.array(z), cmap='RdBu', vmin=np.min(z), vmax=np.max(z))
    # set the limits of the plot to the limits of the data
    cbar = plt.colorbar(c, ax=ax)

    return ax
