import numpy as np, os, sys, pandas as pd

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


def write_lammps_lattice_parameterized(filename="argon_a_0K.in", verbose=False, **kwargs):
    """
    Write LAMMPS input file with customizable parameters
    
    Parameters:
    -----------
    filename : str
        Output filename
    **kwargs : dict
        Parameters to override defaults
    """
    
    # Default parameters
    params = {
        'lattice' : 5.6,
        'tmp_obj' : 1e-5,
        'random_seed': 12345,
        'time_step' : 1e-4,
        'dump_file': 'argon_${lattice}_${tmp_obj}_NPT_*.dat',
        'folder_file' : "${lattice}",
        'num_step' : 20000,
        'press_obj' : 1.01325,
        'npt_thermo_freq' : 0.02,
        'nvt_thermo_freq' : 0.02,
        'num_chains' : 5,
        'num_mtk' : 5,
        'skin' : 2.0,
        'mass' : 39.95,
        'tpotential_lj' : "lj/cut",
    }
    
    # Update with user-provided parameters
    params.update(kwargs)
    
    lammps_content = f"""

    units metal 
    dimension 3
    boundary p p p
    atom_style atomic

    # ----------- Variables ------------
    variable tmp_obj equal {params['tmp_obj']}
    variable time_step equal {params['time_step']}
    variable dump_file string "{params['dump_file']}"
    variable folder_file string "{params['folder_file']}"
    variable num_step equal {params['num_step']}
    variable press_obj equal {params['press_obj']}
    variable npt_thermo_freq equal {params['npt_thermo_freq']}
    variable nvt_thermo_freq equal {params['nvt_thermo_freq']}
    variable num_chains equal {params['num_chains']}
    variable num_mtk equal {params['num_mtk']}
    variable skin equal {params['skin']}
    variable mass equal {params['mass']}
    variable tpotential_lj string "{params['tpotential_lj']}"
    variable lattice equal {params['lattice']}

    # ----------- Lattice and Geometry ------------
    lattice fcc ${{lattice}}
    region box block 0 4 0 4 0 4
    create_box 1 box
    create_atoms 1 box

    # ----------- LJ Potential ------------
    mass 1 ${{mass}}
    velocity all create ${{tmp_obj}} {params['random_seed']} mom yes rot yes dist gaussian
    pair_style ${{tpotential_lj}} 8.5125
    pair_coeff 1 1 0.0103 3.405

    # ----------- Neighbor List ------------
    neighbor ${{skin}} bin
    neigh_modify check yes

    # ----------- Dump Settings ------------
    shell mkdir dump_NVT/${{lattice}}
    compute pe_atom all pe/atom
    compute ke_atom all ke/atom
    variable etot atom c_pe_atom+c_ke_atom
    dump 1 all custom 100 dump_NVT/${{lattice}}/${{dump_file}} id type x y z vx vy vz fx fy fz v_etot

    # ----------- NVT Simulation ------------
    thermo 200
    timestep ${{time_step}}
    fix 1 all nvt temp ${{tmp_obj}} ${{tmp_obj}} ${{nvt_thermo_freq}} tchain ${{num_chains}} tloop ${{num_mtk}}      
    run ${{num_step}}
    unfix 1
    """

    with open(filename, 'w') as f:
        f.write(lammps_content)
    
    if verbose: 
        print(f"LAMMPS input file '{filename}' has been created successfully!")
    return params


