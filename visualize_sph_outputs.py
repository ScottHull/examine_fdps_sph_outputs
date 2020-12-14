import os
from copy import copy
import shutil
import numpy as np
import pandas as pd
import moviepy.editor as mpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.centering import center_of_mass


class BuildMovie:

    def __init__(self, output_path, to_path, num_files, fps=10, start_from=0, new_sim=True, colorize_particles=True,
                 dimension='3', file_name="sph_output.mp4", num_processes=1, focus_process=None, center=True,
                 center_on_target_iron=False):
        self.output_path = output_path
        self.to_path = to_path
        self.file_name = file_name
        self.num_files = num_files
        self.start_file = start_from
        self.curr_file = start_from
        self.fps = fps
        self.colorize_particles = colorize_particles
        self.animation_filename = 'sph_output.mp4'
        self.dimension = dimension
        self.num_processes = num_processes
        self.curr_process = 0
        self.focus_process = focus_process
        self.center = center
        self.__center_on_target_iron = center_on_target_iron
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
        else:
            ax = fig.add_subplot(111)
            x = np.array([])
            y = np.array([])
            z = np.array([])
            particle_id = np.array([])
            colors = np.array([])
            mass = np.array([])
            for proc in range(0, self.num_processes, 1):
                savestate = copy(self.curr_file)
                self.curr_file = file_num
                particle_id_t, x_t, y_t, z_t, colors_t, mass_t = self.__read_sph_file()
                x = np.concatenate((x, x_t))
                y = np.concatenate((y, y_t))
                z = np.concatenate((z, z_t))
                particle_id = np.concatenate((particle_id, particle_id_t))
                colors = np.concatenate((colors, colors_t))
                mass = np.concatenate((mass, mass_t))
                self.curr_file = savestate
                self.curr_process = proc
                print(self.__get_filename())
            if self.center:
                com = center_of_mass(x_coords=x, y_coords=y, z_coords=z, masses=mass, particle_ids=particle_id,
                                     target_iron=self.__center_on_target_iron)
                x = x - com[0]
                y = y - com[1]
                z = z - com[2]
            x = np.array([i for index, i in enumerate(x) if particle_id[index] % 2 == 0] + [i for index, i in enumerate(x) if particle_id[index] % 2 != 0])
            y = np.array([i for index, i in enumerate(y) if particle_id[index] % 2 == 0] + [i for index, i in enumerate(y) if particle_id[index] % 2 != 0])
            z = np.array([i for index, i in enumerate(z) if particle_id[index] % 2 == 0] + [i for index, i in enumerate(z) if particle_id[index] % 2 != 0])
            colors = np.array([i for index, i in enumerate(colors) if particle_id[index] % 2 == 0] + [i for index, i in enumerate(colors) if particle_id[index] % 2 != 0])
            mass = np.array([i for index, i in enumerate(mass) if particle_id[index] % 2 == 0] + [i for index, i in enumerate(mass) if particle_id[index] % 2 != 0])
            particle_id = np.array([i for index, i in enumerate(particle_id) if particle_id[index] % 2 == 0] + [i for index, i in enumerate(particle_id) if particle_id[index] % 2 != 0])
            if self.colorize_particles:
                ax.scatter(x, y, c=colors, alpha=alpha)
            else:
                ax.scatter(x, y, c='black', alpha=alpha)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_xbound(-5e7, 5e7)
            ax.set_ybound(-5e7, 5e7)
            # ax.axis('equal')
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

    def __animate(self, file_name="sph_output.mp4"):
        frames = [self.to_path + "/output_{}.png".format(i) for i in range(self.start_file, self.num_files + 1)]
        animation = mpy.ImageSequenceClip(frames, fps=self.fps, load_images=True)
        animation.write_videofile(file_name, fps=self.fps)

    def build_animation(self, save=False):
        for i in range(self.start_file, self.num_files + 1):
            self.__make_scene(savefig=save, file_num=i)
            self.curr_process = 0
        self.__animate(file_name=self.file_name)
        print("Your animation is now available at {}/{}.".format(os.getcwd(), self.animation_filename))
        return None

    def return_single_scene(self, file_num):
        self.__make_scene(savefig=False, file_num=file_num, alpha=1.0)
        return None

    def build_animation_from_existing(self, path, num_files, start_from_file_number=0):
        num_files_copy = copy(self.num_files)

path = "/Users/scotthull/Desktop/test3"
mov = BuildMovie(
    output_path=path,
    to_path=os.getcwd() + "/sph_visualization",
    num_files=1000,
    fps=30,
    colorize_particles=True,
    dimension='3',
    start_from=0,
    num_processes=20,
    focus_process=None,
    file_name="sph_output.mp4",
    center=True,
    center_on_target_iron=False
)

mov.build_animation(save=True)

# mov.return_single_scene(file_num=0)
