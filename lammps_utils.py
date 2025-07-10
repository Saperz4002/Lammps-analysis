# lammps_utils.py
import numpy as np
import os

def extract_data(filename):
    with open(filename, encoding='ISO-8859-1') as file:
        data = file.readlines()

    time_step = None
    x, y, z = [], [], []
    vx, vy, vz = [], [], []
    fx, fy, fz = [], [], []
    etot = []
    atom_data_started = False

    for i, line in enumerate(data):
        if line.startswith("ITEM: TIMESTEP"):
            time_step = int(data[i + 1].strip())
        elif line.startswith("ITEM: ATOMS"):
            atom_data_started = True
            continue
        if atom_data_started:
            if line.strip() == "":
                break
            components = list(map(float, line.split()[2:]))
            x.append(components[0])
            y.append(components[1])
            z.append(components[2])
            vx.append(components[3])
            vy.append(components[4])
            vz.append(components[5])
            fx.append(components[6])
            fy.append(components[7])
            fz.append(components[8])
            etot.append(components[9])

    return time_step, x, y, z, vx, vy, vz, fx, fy, fz, etot

def etotal(e_vec):
    return sum(e_vec)

def get_dump_files(directory, prefix="argon.lj."):
    return [os.path.join(directory, f) for f in os.listdir(directory) if f.startswith(prefix)]

def load_all_data(dump_files):
    results = []
    all_velocities = []

    for filepath in dump_files:
        time_step, x, y, z, vx, vy, vz, fx, fy, fz, etot = extract_data(filepath)
        if etotal(etot) == 0:
            continue
        total_energy = etotal(etot)
        results.append((time_step, total_energy))
        all_velocities.append(np.array(vx))
    
    return results, all_velocities
