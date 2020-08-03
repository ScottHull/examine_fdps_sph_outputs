import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np

def plot_heatmap(x, y, z, a, b, center):
    ax = plt.figure().add_subplot(111)
    sc = ax.scatter(x, y, c=z, marker="+")
    e = Ellipse(xy=(0, 0), width=a * 2.0, height=b * 2.0, alpha=0.2, color="blue")
    ax.scatter(center[0], center[1], marker="*", color="red", s=60)
    ax.add_artist(e)
    cbar = plt.colorbar(sc, ax=ax)

    return ax
