import matplotlib.pyplot as plt

def plot_heatmap(x, y, z):
    ax = plt.figure().add_subplot(111)
    c = ax.pcolormesh([x, y], z, cmap='RdBu', vmin=min(z), vmax=max(z))
    # set the limits of the plot to the limits of the data
    cbar = plt.colorbar(c, ax=ax)

    return ax
