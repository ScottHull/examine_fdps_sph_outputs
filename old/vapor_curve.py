import matplotlib.pyplot as plt
import pandas as pd


class PhaseSpace:

    def __init__(self, phase="duniteS"):
        self.phase = phase
        self.phase_df = pd.DataFrame({})

        if phase == "duniteS":
            self.phase_df = pd.read_fwf("src/phase_data/duniteS_vapour_curve.txt", skiprows=1,
                                        names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                               "entropy_sol_liq", "entropy_vap"])
        elif phase == "dunite4":
            self.phase_df = pd.read_fwf("src/phase_data/duniteS_vapour_curve.txt", skiprows=1,
                                        names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                               "entropy_sol_liq", "entropy_vap"])

    def plot_phase_space(self, x_header="entropy", y_header="temperature"):
        ax = plt.figure().add_subplot(111)
        if x_header == "entropy":
            entropies = ['entropy_sol_liq', 'entropy_vap']
            entropies_labels = ["Solid-Liquid Entropy Phase Boundary", "Liquid-Vapor Entropy Phase Boundary"]
            for index, i in enumerate(entropies):
                label = entropies_labels[index]
                ax.plot(self.phase_df[i], self.phase_df[y_header], linewidth=2.0, label=i)
            ax.fill_betweenx(self.phase_df[y_header], self.phase_df[entropies[0]], self.phase_df[entropies[1]],
                             alpha=0.6, color='aqua')
        else:
            ax.plot(self.phase_df[x_header], self.phase_df[y_header], linewidth=2.0, label=x_header)
        ax.set_xlabel(x_header)
        ax.set_ylabel(y_header)
        ax.set_title(self.phase)
        ax.grid()
        ax.legend()
        plt.show()


PhaseSpace().plot_phase_space()
