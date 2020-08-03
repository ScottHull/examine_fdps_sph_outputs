import matplotlib.pyplot as plt
import numpy as np

def plot_heatmap(x, y, z):
    ax = plt.figure().add_subplot(111)
    m = np.meshgrid(x, y)
    c = ax.pcolormesh(m, z, cmap='RdBu', vmin=min(z), vmax=max(z))
    # set the limits of the plot to the limits of the data
    cbar = plt.colorbar(c, ax=ax)

    return ax
