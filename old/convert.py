import glob
import math
import sys

G = 6.67e-11
R = 6400e+3
M = 6.0e+24

for s in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    print(s)
    f = open("{0:0>5}.pov".format(s), 'w')
    f.write('#include "colors.inc"\n')
    f.write('#include "stones.inc"\n')
    f.write('#include "metals.inc"\n')
    f.write('#include "textures.inc"\n')
    f.write('#include "glass.inc"\n')
    f.write('#include "glass_old.inc"\n')
    f.write('#include "stars.inc"\n')

    f.write('camera{\n')
    f.write('\tlocation <0., 0., -20.0>\n')
    f.write('\tlook_at <0, 0, 0>\n')
    f.write('\tright x * image_width/image_height\n')
    f.write('}\n')

    f.write('light_source{\n')
    f.write('\t<0.0, 0.0, -20.0>\n')
    f.write('\tcolor rgb 1.\n')
    f.write('\tshadowless\n')
    f.write('}\n')

    # f.write('object{')
    # f.write('\tsphere{')
    # f.write('\t\t<0,0,0>, 1')
    # f.write('\t\tinverse')
    # f.write('\t\ttexture{')
    # f.write('\t\t\tStarfield1')
    # f.write('\t\t}')
    # f.write('\tscale 1000')
    # f.write('\t}')
    # f.write('}')

    for i in range(32):
        with open("result/{0:0>5}_00032_{1:0>5}.dat".format(s, i), 'r') as file:
            file.readline()
            file.readline()
            for line in file:
                data = line.split('\t')

                id = int(data[0])
                mat = int(data[1])
                dens = float(data[9])
                mass = float(data[2])
                radius = 1.2 * math.pow(mass / dens, 1. / 3.) / R
                rx = float(data[3]) / R
                ry = float(data[4]) / R
                rz = float(data[5]) / R
                if rz < 0:
                    continue
                f.write('object{\n')

                f.write('\tsphere{\n')
                f.write('\t\t<{x}, {y}, {z}>, {r}\n'.format(x=rx, y=ry, z=rz, r=radius))
                f.writfrom
                mpl_toolkits.mplot3d
                import Axes3De

                ('\t}\n')
                if mat == 0:
                    texture = "T_Chrome_1A"
                    f.write('\ttexture{\n')
                    f.write('\t\t' + texture + '\n')
                    f.write('\t}\n')
                elif mat == 1:
                    '''
                    texture = "T_Stone28"
                    f.write('\ttexture{\n')
                    f.write('\t\t' + texture + '\n')
                    f.write('\t}\n')
                    '''
                    texture = "M_Orange_Glass"
                    f.write('\tmaterial{\n')
                    f.write('\t\t' + texture + '\n')
                    f.write('\t}\n')
                    '''
                    f.write('\tmaterial{\n')
                    f.write('\t\ttexture{\n')
                    f.write('\t\t\tpigment{\n')
                    f.write('\t\t\t\tcolor VeryDarkBrown\n')
                    f.write('\t\t\t}\n')
                    f.write('\t\t\tfinish{\n')
                    f.write('\t\t\t\tF_Glass8\n')
                    f.write('\t\t\t}\n')
                    f.write('\t\t}\n')
                    f.write('\t\tinterior{\n')
                    f.write('\t\t\tI_Glass_Caustics2\n')
                    f.write('\t\t\tior 0\n')
                    f.write('\t\t\tfade_color Col_Fluorite_04\n')
                    f.write('\t\t}\n')
                    f.write('\t}\n')
                    '''

                elif mat == 2:
                    '''
                    texture = "M_Ruby_Glass"
                    f.write('\tmaterial{\n')
                    f.write('\t\t' + texture + '\n')
                    f.write('\t}\n')
                    '''
                    f.write('\tmaterial{\n')
                    f.write('\t\ttexture{\n')
                    f.write('\t\t\tpigment{\n')
                    f.write('\t\t\t\tcolor Col_Glass_Ruby\n')
                    f.write('\t\t\t}\n')
                    f.write('\t\t\tfinish{\n')
                    f.write('\t\t\t\tF_Glass2\n')
                    f.write('\t\t\t}\n')
                    f.write('\t\t}\n')
                    f.write('\t\tinterior{\n')
                    f.write('\t\t\tI_Glass_Caustics1\n')
                    f.write('\t\t\tior 0\n')
                    f.write('\t\t\tfade_color Col_Fluorite_04\n')
                    f.write('\t\t}\n')
                    f.write('\t}\n')
                elif mat == 3:
                    texture = "Silver_Texture"
                    f.write('\ttexture{\n')
                    # f.write('\t\tpigment{\n')
                    # f.write('\t\t\tcolor DimGrey\n')
                    # f.write('\t\t}\n')
                    f.write('\t\t' + texture + '\n')
                    f.write('\t}\n')
                elif mat == 4:
                    '''
                    texture = "T_Stone41"
                    f.write('\ttexture{\n')
                    f.write('\t\t' + texture + '\n')
                    f.write('\t}\n')
                    '''
                    texture = "M_Dark_Green_Glass"
                    f.write('\ttexture{\n')
                    f.write('\t\tpigment{\n')
                    f.write('\t\t\tcolor SlateBlue\n')
                    f.write('\t\t}\n')
                    # f.write('\t\t' + texture + '\n')
                    f.write('\t}\n')
                f.write('\tfinish{\n')
                f.write('\t\tambient .5\n')
                f.write('\t\tdiffuse 1.0\n')
                f.write('\t\tphong 0\n')
                f.write('\t}\n')
                f.write('\tnormal{\n')
                f.write('\t\tagate 0.8\n')
                f.write('\t\tscale 0.25\n')
                f.write('\t}\n')
                f.write('}\n')
