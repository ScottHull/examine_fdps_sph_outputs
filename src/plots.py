import matplotlib.pyplot as plt
import numpy as np

def plot_heatmap(x, y, z):
    ax = plt.figure().add_subplot(111)
    # y, x = np.meshgrid(x, y)
    x = np.unique(x)
    y = np.unique(y)
    X, Y = np.meshgrid(x, y)
    Z = np.array(z).reshape(len(y), len(x))
    c = ax.pcolormesh(X, Y, Z)
    cbar = plt.colorbar(c, ax=ax)

    return ax
