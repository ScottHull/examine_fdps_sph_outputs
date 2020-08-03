import pandas as pd
import matplotlib.pyplot as plt

iron___ = pd.read_fwf("/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/FDPS_SPH/eos/iron___.rho_u.txt",
                      skiprows=2, header=None)
ironS2 = pd.read_fwf("/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/FDPS_SPH/eos/ironS2.rho_u.txt",
                     skiprows=2, header=None)
ironC = pd.read_fwf("/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/FDPS_SPH/eos/ironC.rho_u.txt",
                     skiprows=2, header=None)

ax1 = plt.figure().add_subplot(111)
ax1.scatter(iron___[5], iron___[1], marker="+", color='blue', label="iron___")
# ax1.scatter(ironS2[5], ironS2[1], marker="+", color='red', label="ironS2")
ax1.scatter(ironC[5], ironC[1], marker="+", color='green', label="ironC")
ax1.set_xlabel("Entropy")
ax1.set_ylabel("Internal Energy")
ax1.set_title("Entropy vs. Internal Energy")
ax1.grid()
ax1.legend()

ax2 = plt.figure().add_subplot(111)
ax2.scatter(iron___[5], iron___[1], marker="+", color='blue', label="iron___")
# ax2.scatter(ironS2[5], ironS2[1], marker="+", color='red', label="ironS2")
ax2.scatter(ironC[5], ironC[1], marker="+", color='green', label="ironC")
ax2.set_xlabel("Entropy")
ax2.set_ylabel("Temperature")
ax2.set_title("Entropy vs. Temperature")
ax2.grid()
ax2.legend()

ax3 = plt.figure().add_subplot(111)
ax3.scatter(iron___[5], iron___[3], marker="+", color='blue', label="iron___")
# ax3.scatter(ironS2[5], ironS2[3], marker="+", color='red', label="ironS2")
ax3.scatter(ironC[5], ironC[3], marker="+", color='green', label="ironC")
ax3.set_xlabel("Entropy")
ax3.set_ylabel("Pressure")
ax3.set_title("Entropy vs. Pressure")
ax3.grid()
ax3.legend()

ax4 = plt.figure().add_subplot(111)
ax4.scatter(iron___[5], iron___[4], marker="+", color='blue', label="iron___")
# ax4.scatter(ironS2[5], ironS2[4], marker="+", color='red', label="ironS2")
ax4.scatter(ironC[5], ironC[4], marker="+", color='green', label="ironC")
ax4.set_xlabel("Entropy")
ax4.set_ylabel("Sound Speed")
ax4.set_title("Entropy vs. Sound Speed")
ax4.grid()
ax4.legend()

ax5 = plt.figure().add_subplot(111)
ax5.scatter(iron___[5], iron___[5], marker="+", color='blue', label="iron___")
# ax5.scatter(ironS2[5], ironS2[5], marker="+", color='red', label="ironS2")
ax5.scatter(ironC[5], ironC[5], marker="+", color='green', label="ironC")
ax5.set_xlabel("Entropy")
ax5.set_ylabel("Entropy")
ax5.set_title("Entropy vs. Entropy")
ax5.grid()
ax5.legend()

plt.show()
