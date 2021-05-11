from src.interpolation import GenericTrilinearInterpolation
from scipy.integrate import odeint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def __adiabat_eq(T, y, c_p, alpha_v, g):
    # dT/dy = (alpha_v g T) / c_p
    dTdy = (alpha_v * g * T) / c_p
    return dTdy


def __hydrostatic_eq(rho, g):
    dPdy = rho * g
    return dPdy


def adiabat(T_0, c_p, g, alpha_v, y):
    return [i[0] for i in odeint(__adiabat_eq, T_0, y, args=(c_p, alpha_v, g))]


def hydrostatic(P_0, rho, g, y):
    return [i[0] for i in odeint(__hydrostatic_eq, P_0, y, args=(rho, g))]


def map_hydrostatic_pressure_to_entropy(hydrostatic_pressures):
    pass


def manual_adiabat(T_0, c_p, g, alpha_v, y_resolution, max_y):
    T = T_0
    y = 0
    temps = [T_0]
    depths = [y]
    while y <= max_y:
        y += y_resolution
        dTdy = (alpha_v * g * T) / c_p
        dT = dTdy * y_resolution
        T += dT
        temps.append(T)
        depths.append(y)
    return temps, depths


min_depth = 0
max_depth = 7000 * 1000
depth_resolution = int((max_depth / 1000) + 1)
depth_interval = np.linspace(min_depth, max_depth + depth_resolution, depth_resolution)
T_0 = 1600
c_p = 1000
alpha_v = 3 * 10 ** -5
g = 9.8
a = adiabat(T_0=T_0, c_p=c_p, alpha_v=alpha_v, g=g, y=depth_interval)
manual_a, manual_y = manual_adiabat(T_0=T_0, c_p=c_p, alpha_v=alpha_v, g=g, max_y=max_depth,
                                    y_resolution=depth_resolution)

ax = plt.figure().add_subplot(111)
ax.plot(list(a), [i / 1000 for i in depth_interval], color='blue', linewidth=2.0, label="ODE Solver")
ax.plot(manual_a, [i / 1000 for i in manual_y], color='red', linewidth=2.0, linestyle="--", label="Manual")
ax.invert_yaxis()
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Depth (km)")
ax.set_title("Adidbatic Gradient (max depth: {} km)".format(max_depth / 1000))
ax.grid()
ax.legend()

plt.show()
