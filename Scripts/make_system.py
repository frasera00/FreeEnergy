#setup
from Building import Builder
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import pyplot
import MDAnalysis as mda
import tqdm as tqdm
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import sys
import os


def write_gromacs_top(system_name, molecules, output_file, includes):
    with open(output_file, 'w') as f:
        for include in includes:
            f.write(f'#include "{include}"\n')
        
        f.write(f'[ system ]\n{system_name}\n\n')
        
        f.write('[ molecules ]\n')
        for molecule, count in molecules.items():
            f.write(f'  {molecule}     {count}\n')

def insert_gas_with_radius(file, insert_file, outfile, R, resname, N, outplot=None):
    u = mda.Universe(file)
    bulk = u.select_atoms("resname BULK")
    pis = u.select_atoms("resname PIS")
    pis1 = u.select_atoms("resname PIS1")

    # Define the bounds for molecule placement
    x_min, x_max = pis.positions[:,0].min(), pis.positions[:,0].max()
    y_min, y_max = pis.positions[:,1].min(), pis.positions[:,1].max()
    x_min, x_max = x_min + 2 * R, x_max - 2 * R
    y_min, y_max = y_min + 2 * R, y_max - 2 * R

    N = N // 2  # Half molecules per region

    inserted_centroids = [] 
    atoms_to_replace = []
    atoms_to_insert = []

    fig, ax = plt.subplots()

    pbar = tqdm.tqdm(total=2*N, desc="Inserting molecules")
    for mins in [(pis.positions[:,2].max(), bulk.positions[:,2].min()), 
                 (bulk.positions[:,2].max(), pis1.positions[:,2].min())]:
        z_min, z_max = mins
        z_min, z_max = z_min + 2 * R, z_max - 2 * R
        water_atoms = u.select_atoms(f"(resname W or resname SW or resname TW) and prop z >= {z_min} and prop z <= {z_max} and prop x >= {x_min} and prop x <= {x_max} and prop y >= {y_min} and prop y <= {y_max}")

        for _ in range(N):
            pbar.update(1)
            max_attempts = 1000
            for _ in range(max_attempts):
                centroid = np.array([
                    np.random.uniform(x_min, x_max),
                    np.random.uniform(y_min, y_max),
                    np.random.uniform(z_min, z_max)
                ])
                insert_u = mda.Universe(insert_file)

                if all(np.linalg.norm(centroid - c) >= 2 * R for c in inserted_centroids):
                    inserted_centroids.append(centroid)
                    
                    distances = np.linalg.norm(water_atoms.positions - centroid, axis=1)
                    atoms_to_replace.append(water_atoms[distances <= R])

                    insert_centroid = insert_u.atoms.center_of_geometry()
                    shift = centroid - insert_centroid
                    insert_u.atoms.translate(shift)  # Apply translation correctly

                    atoms_to_insert.append(insert_u.atoms)
                   
                    #print(f"Inserted {resname} at {centroid}, replaced {len(atoms_to_replace)} water molecules.")
                    break
            else:
                print(f"Failed to place molecule {len(inserted_centroids) + 1} without overlap after {max_attempts} attempts.")

    remaining_atoms = u.atoms
    for group in atoms_to_replace:
        remaining_atoms = remaining_atoms.difference(group)

    all_atoms = [remaining_atoms]
    for group in atoms_to_insert:
        all_atoms.append(group)

    new_universe = mda.Merge(*all_atoms)
    new_universe.dimensions = u.dimensions
    new_universe.atoms.write(outfile)

    water = new_universe.select_atoms("resname W or resname SW or resname TW")
    gas = new_universe.select_atoms(f"resname {resname}")

    ax.scatter(water.positions[:, 2], water.positions[:, 0], s=0.1, alpha=0.5, label="Water")
    ax.scatter(gas.positions[:, 2], gas.positions[:, 0], color='red', s=50, label=f"{resname}")
    ax.legend()

    if outplot:
        fig.savefig(outplot, dpi=200)
    else:
        plt.clf()


root = '/scratch/rasera/FreeEnergy'
mol = sys.argv[1]
resname = sys.argv[2]
ngas = int(sys.argv[3])
outfile = f"{root}/start_files/eqruns/{mol}/{ngas}/FE_{resname}_{ngas}.gro"
outplot = f"{root}/start_files/eqruns/{mol}/{ngas}/FE_{resname}_{ngas}.png"
insert_file = f"{root}/start_files/{mol}.gro"
file = f"{root}/start_files/FE.gro"
R = 5 # Å

# Create a single new folder
folder_name = root+f"/start_files/eqruns/{mol}/{ngas}"
os.makedirs(folder_name, exist_ok=True) 

print(f"Insert {ngas} molecules, res {resname}, in {mol} system")
insert_gas_with_radius(file, insert_file, outfile, R, resname, N=ngas)

u = mda.Universe(root+f"/start_files/eqruns/{mol}/{ngas}/FE_{resname}_{ngas}.gro")
SIL = u.select_atoms("name N1L1").n_atoms
LAY = u.select_atoms("name N1L").n_atoms
PIS = u.select_atoms("resname PIS").n_atoms
PIS1 = u.select_atoms("resname PIS1").n_atoms
BULK = u.select_atoms("name N1B").n_atoms
GASMOL = u.select_atoms(f"resname {resname}").n_atoms
W = u.select_atoms("resname W").n_atoms
SW = u.select_atoms("resname SW").n_atoms
TW = u.select_atoms("resname TW").n_atoms

print(SIL,LAY,PIS,PIS1,BULK,GASMOL)

molecules = {
    "PIS": PIS,
    "SIL": SIL,
    "LAY": LAY,
    "BULK": BULK,
    "PIS1": PIS,
    "TW": TW,
    "SW": SW,
    "W": W,
    f"{resname}": GASMOL
}

includes = [
    "martini_v3.0.0_N1Lay.itp",
    "martini_v3.0.0_solvents_v1.itp",
    "lay1-sil.itp",
    "surf.itp",
]

write_gromacs_top(system_name="System",
                  molecules=molecules,
                  output_file=root+f"/start_files/eqruns/{mol}/{ngas}/FE_{resname}_{ngas}.top",
                  includes=includes)