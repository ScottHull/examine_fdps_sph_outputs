import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats.kde import gaussian_kde
import numpy as np

def plot_heatmap(x, y, z, a, b, center):
    ax = plt.figure().add_subplot(111)
    sc = ax.scatter(x, y, c=z, marker="+")
    e = Ellipse(xy=(0, 0), width=a * 2.0, height=b * 2.0, alpha=0.2, color="blue")
    ax.scatter(center[0], center[1], marker="*", color="red", s=200)
    ax.add_artist(e)
    cbar = plt.colorbar(sc, ax=ax)

    return ax

def plot_particle_density_heatmap(x, y):
    x = np.array(x)
    y = np.array(y)
    k = gaussian_kde(np.vstack([x, y]))
    xi, yi = np.mgrid[x.min():x.max():x.size ** 0.5 * 1j, y.min():y.max():y.size ** 0.5 * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    ax = plt.figure().add_subplot(111)
    # alpha=0.5 will make the plots semitransparent
    ax.pcolormesh(xi, yi, zi.reshape(xi.shape), alpha=0.5)

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())

    return ax
