# main.py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
from lammps_utils import get_dump_files, load_all_data

if __name__ == "__main__":
    directory = os.getcwd()
    dump_files = get_dump_files(directory, prefix="argon.lj.")
    results, all_velocities = load_all_data(dump_files)

    time_steps = [result[0] for result in results]

    ####Histogram of velocities
    hist_matrix = []
    all_v_flat = np.concatenate(all_velocities)
    vmin = np.min(all_v_flat)
    vmax = np.max(all_v_flat)
    step_interval = 100
    for v in all_velocities[::step_interval]:
        counts, _ = np.histogram(v, bins=20, range=(vmin, vmax))
        hist_matrix.append(counts)
    hist_matrix = np.array(hist_matrix)

    fig = plt.figure(figsize=(12,7))
    ax = fig.add_subplot(111, projection='3d')
    T, B = hist_matrix.shape
    xpos, ypos = np.meshgrid(np.arange(B), np.arange(T))
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros_like(xpos)

    dx = dy = 0.8
    dz = hist_matrix.flatten()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, shade=True)
    ax.set_xlabel('Velocity')
    ax.set_ylabel('Time Step / 100')
    ax.set_zlabel('Frequency')
    ax.set_title('Velocities Over Time')

    plt.show()
