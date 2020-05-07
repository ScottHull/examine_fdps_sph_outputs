import os
from copy import copy
import shutil
import pandas as pd
import moviepy.editor as mpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class BuildMovie:

    def __init__(self, output_path, to_path, num_files, fps=10, start_from=0, new_sim=True, colorize_particles=True,
                 dimension='3'):
        self.output_path = output_path
        self.to_path = to_path
        self.num_files = num_files
        self.start_file = start_from
        self.curr_file = start_from
        self.fps = fps
        self.colorize_particles = colorize_particles
        self.animation_filename = 'sph_output.mp4'
        self.file_format = "results.{}_00001_00000.dat"
        self.dimension = dimension
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
        return self.output_path + "/" + self.file_format.format(str(self.curr_file).zfill(5))

    def __read_sph_file(self):
        df = pd.read_csv(self.__get_filename(), sep='\t', skiprows=2, header=None)
        particle_id = df[1]
        x = df[3]
        y = df[4]
        z = df[5]
        colors = []
        if self.colorize_particles:
            colors = [self.color_map[int(i)] for i in particle_id]
        return particle_id, x, y, z, colors

    def __make_scene(self, savefig=True, file_num=None, alpha=1.0):
        if savefig:
            particle_id, x, y, z, colors = self.__read_sph_file()
        else:
            savestate = copy(self.curr_file)
            self.curr_file = file_num
            particle_id, x, y, z, colors = self.__read_sph_file()
            self.curr_file = savestate
        fig = plt.figure()
        fig.set_size_inches(18.5, 10.5)
        if self.dimension == '3':
            ax = Axes3D(fig)
            if self.colorize_particles:
                ax.scatter(x, y, z, c=colors, alpha=alpha)
            else:
                ax.scatter(x, y, z, c='black', alpha=alpha)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            # ax.set_xbound(-5e7, 5e7)
            # ax.set_ybound(-5e7, 5e7)
            # ax.set_zbound(-5e7, 5e7)
        else:
            ax = fig.add_subplot(111)
            if self.colorize_particles:
                ax.scatter(x, y, c=colors, alpha=alpha)
            else:
                ax.scatter(x, y, c='black', alpha=alpha)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            # ax.set_xbound(-5e7, 5e7)
            # ax.set_ybound(-5e7, 5e7)
            ax.axis('equal')
        if savefig:
            ax.set_title("Iteration: {}".format(self.curr_file))
            print("Built scene: {}".format(self.curr_file))
            fig.savefig(self.to_path + "/output_{}.png".format(self.curr_file), format='png', dpi=100)
            self.curr_file += 1
        else:
            plt.show()

    def __animate(self):
        frames = [self.to_path + "/output_{}.png".format(i) for i in range(self.start_file, self.num_files + 1)]
        animation = mpy.ImageSequenceClip(frames, fps=self.fps, load_images=True)
        animation.write_videofile('sph_output.mp4', fps=self.fps)

    def build_animation(self):
        for i in range(self.start_file, self.num_files + 1):
            self.__make_scene(savefig=True, file_num=None)
        self.__animate()
        print("Your animation is now available at {}/{}.".format(os.getcwd(), self.animation_filename))
        return None

    def return_single_scene(self, file_num):
        self.__make_scene(savefig=False, file_num=file_num, alpha=1.0)
        return None

    def build_animation_from_existing(self, path, num_files, start_from_file_number=0):
        num_files_copy = copy(self.num_files)

path = "/Users/scotthull/Desktop/test2"
mov = BuildMovie(
    output_path=path,
    to_path=os.getcwd() + "/sph_visualization",
    num_files=1000,
    fps=30,
    colorize_particles=True,
    dimension='3',
    start_from=0
)

mov.build_animation()

# mov.return_single_scene(file_num=0)
