import matplotlib.pyplot as plt
import numpy as np

def plot_heatmap(x, y, z):
    ax = plt.figure().add_subplot(111)
    sc = ax.scatter(x, y, c=z, marker="+")
    cbar = plt.colorbar(sc, ax=ax)

    return ax
