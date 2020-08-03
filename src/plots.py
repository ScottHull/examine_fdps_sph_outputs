import matplotlib.pyplot as plt
import numpy as np

def plot_heatmap(x, y, z):
    ax = plt.tricontour(x, y, z, 15, linewidths=0.5, colors='k')
    return ax
