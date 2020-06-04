import os
from operator import add
from copy import copy
import shutil
import pandas as pd
import numpy as np
import moviepy.editor as mpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class BuildMovie:

    def __init__(self, output_path, to_path, num_files, fps=10, start_from=0, new_sim=True, colorize_particles=True,
                 dimension='3', num_processes=1, focus_process=None, center=False):
        self.output_path = output_path
        self.to_path = to_path
        self.num_files = num_files
        self.start_file = start_from
        self.curr_file = start_from
        self.fps = fps
        self.center = center
        self.colorize_particles = colorize_particles
        self.animation_filename = 'sph_output.mp4'
        self.dimension = dimension
        self.num_processes = num_processes
        self.curr_process = 0
        self.focus_process = focus_process
        self.file_format = "results.{}_{}_{}.dat"
        self.color_map = {
            0: 'red',
            1: 'blue',
            2: 'yellow',
            3: 'purple',
        }
        if new_sim is True:
            if os.path.exists(self.to_path):
                shutil.rmtree(self.to_path)
            os.mkdir(self.to_path)

    def __get_filename(self):
        return self.output_path + "/" + self.file_format.format(str(self.curr_file).zfill(5),
                                                                str(self.num_processes).zfill(5),
                                                                str(self.curr_process).zfill(5))

    def __read_sph_file(self):
        df = pd.read_csv(self.__get_filename(), sep='\t', skiprows=2, header=None)
        particle_id = df[1]
        mass = df[2]
        x = df[3]
        y = df[4]
        z = df[5]
        colors = []
        if self.colorize_particles:
            colors = [self.color_map[int(i)] for i in particle_id]
        return particle_id, x, y, z, colors, mass

    def __make_scene(self, savefig=True, file_num=None, alpha=1.0):
        fig = plt.figure()
        fig.set_size_inches(18.5, 10.5)
        if self.dimension == '3':
            ax = Axes3D(fig)
            if self.focus_process is None:
                for proc in range(0, self.num_processes, 1):
                    self.curr_process = proc
                    if savefig:
                        particle_id, x, y, z, colors, mass = self.__read_sph_file()
                    else:
                        savestate = copy(self.curr_file)
                        self.curr_file = file_num
                        particle_id, x, y, z, colors, mass = self.__read_sph_file()
                        self.curr_file = savestate
                    print(self.__get_filename())
                    if self.colorize_particles:
                        ax.scatter(x, y, z, c=colors, alpha=alpha)
                    else:
                        ax.scatter(x, y, z, c='black', alpha=alpha)
                    self.curr_process = proc
            else:
                self.curr_process = self.focus_process
                if savefig:
                    particle_id, x, y, z, colors, mass = self.__read_sph_file()
                else:
                    savestate = copy(self.curr_file)
                    self.curr_file = file_num
                    particle_id, x, y, z, colors, mass = self.__read_sph_file()
                    self.curr_file = savestate
                if self.colorize_particles:
                    ax.scatter(x, y, z, c=colors, alpha=alpha)
                else:
                    ax.scatter(x, y, z, c='black', alpha=alpha)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            ax.set_xbound(-9e6, 9e6)
            ax.set_ybound(-9e6, 9e6)
            ax.set_zbound(-9e6, 9e6)
            if savefig:
                ax.set_title("Iteration: {}".format(self.curr_file))
                print("Built scene: {}".format(self.curr_file))
                fig.savefig(self.to_path + "/output_{}.png".format(self.curr_file), format='png', dpi=100)
                self.curr_file += 1
            else:
                plt.show()
        else:
            ax = fig.add_subplot(111)
            x_array = np.array([])
            y_array = np.array([])
            mass_array = np.array([])
            colors_array = np.array([])
            for proc in range(0, self.num_processes, 1):
                savestate = copy(self.curr_file)
                self.curr_file = savestate
                self.curr_process = proc
                self.curr_file = file_num
                particle_id, x, y, z, colors, mass = self.__read_sph_file()
                x_array = np.concatenate([x_array, x])
                y_array = np.concatenate([y_array, y])
                mass_array = np.concatenate([mass_array, mass])
                colors_array = np.concatenate([colors_array, colors])
                print(self.__get_filename())
            if self.center:
                x_center, y_center, z_center = self.center_of_mass(x_coords=x_array, y_coords=y_array,
                                                                   z_coords=np.array([]), masses=mass_array)
                x_array -= x_center
                y_array -= y_center
            if self.colorize_particles:
                ax.scatter(x_array, y_array, c=colors_array, alpha=alpha)
            else:
                ax.scatter(x_array, y_array, c='black', alpha=alpha)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            # if self.center:
            #     x_center, y_center, z_center = self.center_of_mass(x_coords=x_array, y_coords=y_array, z_coords=np.array([]), masses=mass_array)
            #     ax.set_xbound(x_center - 5e7, x_center + 5e7)
            #     ax.set_ybound(y_center - 5e7, y_center + 5e7)
            ax.axis('equal')
            if savefig:
                ax.set_title("Iteration: {}".format(self.curr_file))
                print("Built scene: {}".format(self.curr_file))
                fig.savefig(self.to_path + "/output_{}.png".format(self.curr_file), format='png', dpi=100)
                self.curr_file += 1
            else:
                plt.show()

        fig.clear()
        plt.close(fig)

        self.curr_process = 0

    def center_of_mass(self, x_coords, y_coords, z_coords, masses):
        masses = np.array(masses, dtype=np.float32)
        total_mass = float(np.sum(masses))
        x_center = sum([a * b for a, b in zip(x_coords, masses)]) / total_mass
        y_center = sum([a * b for a, b in zip(y_coords, masses)]) / total_mass
        z_center = sum([a * b for a, b in zip(z_coords, masses)]) / total_mass
        return (x_center, y_center, z_center)

    def __animate(self):
        frames = [self.to_path + "/output_{}.png".format(i) for i in range(self.start_file, self.num_files + 1)]
        animation = mpy.ImageSequenceClip(frames, fps=self.fps, load_images=True)
        animation.write_videofile('sph_output.mp4', fps=self.fps)

    def build_animation(self, save=False):
        for i in range(self.start_file, self.num_files + 1):
            self.__make_scene(savefig=save, file_num=i)
            self.curr_process = 0
        self.__animate()
        print("Your animation is now available at {}/{}.".format(os.getcwd(), self.animation_filename))
        return None

    def return_single_scene(self, file_num):
        self.__make_scene(savefig=True, file_num=file_num, alpha=1.0)
        return None

    def build_animation_from_existing(self, path, num_files, start_from_file_number=0):
        num_files_copy = copy(self.num_files)

path = "/Users/scotthull/Desktop/GI_small"
mov = BuildMovie(
    output_path=path,
    to_path=os.getcwd() + "/sph_visualization",
    num_files=4,
    fps=30,
    colorize_particles=True,
    dimension='2',
    start_from=0,
    num_processes=4,
    focus_process=None,
    center=True
)

mov.build_animation(save=True)

# mov.return_single_scene(file_num=0)
