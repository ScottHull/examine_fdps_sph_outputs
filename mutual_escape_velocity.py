from math import sqrt
import pandas as pd

target_path = "/Users/scotthull/Desktop/tar.dat"
impactor_path = "/Users/scotthull/Desktop/imp.dat"

tar_df = pd.read_csv(target_path, skiprows=2, header=None, delimiter="\t")
imp_df = pd.read_csv(impactor_path, skiprows=2, header=None, delimiter="\t")

def calculate_mutual_escape_velocity(m_tar, m_imp, r_tar, r_imp):
    G = 6.674 * 10 ** -11
    return sqrt(2 * G * ((m_tar + m_imp) / (r_tar + r_imp)))

max_x_tar = max([float(i) for i in tar_df[3]])
max_y_tar = max([float(i) for i in tar_df[4]])
max_z_tar = max([float(i) for i in tar_df[5]])
max_x_imp = max([float(i) for i in imp_df[3]])
max_y_imp = max([float(i) for i in imp_df[4]])
max_z_imp = max([float(i) for i in imp_df[5]])
tar_mass = sum([float(i) for i in tar_df[2]])
imp_mass = sum([float(i) for i in imp_df[2]])

v_esc_mutual = calculate_mutual_escape_velocity(m_tar=tar_mass, m_imp=imp_mass, r_imp=max_x_imp, r_tar=max_x_tar)

print(
    "MAX X TAR: {}\n"
    "MAX Y TAR: {}\n"
    "MAX Z TAR: {}\n"
    "MAX X IMP: {}\n"
    "MAX Y IMP: {}\n"
    "MAX Z IMP: {}\n"
    "TAR MASS: {}\n"
    "IMP MASS: {}\n"
    "MUTUAL ESC VELOCITY: {}".format(
        max_x_tar, max_y_tar, max_z_tar,
        max_x_imp, max_y_imp, max_z_imp,
        tar_mass, imp_mass,
        v_esc_mutual
    )
)

print(list(tar_df[2])[0])
print(list(imp_df[2])[0])

