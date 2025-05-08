import numpy as np
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm
import glob

ROOT = "/scratch/rasera"
GASES = {"Chloroform":"CLF",
         "Diethyl_ether":"ETH",
         "Ethylene":"ETY",
         "No_gas":"No_gas",
         "Nitrogen":"NI2",
         "Fullerene":"F15",
         "Xenon":"XEN",
         "Oxygen":"OXY",
         "Bare":"Bare"}


import sys as sys

SYS = "C8,1.2"
INC = "_inc_50"
KAPPA = 0.01
MOL = sys.argv[1]
NGAS = sys.argv[2]

REST_TIME = 50000; EQ_TIME = 500000
N_STEPS = REST_TIME + EQ_TIME; END_STEPS = 550000
AT0 = 8200; ATF = 10800; INCREMENT_AT = 50
TARGETS = np.arange(AT0, ATF+INCREMENT_AT, INCREMENT_AT)

GAS = GASES[MOL]
data = []

for i,TARGET in enumerate(tqdm(TARGETS)):
    try:
        gro = glob.glob(f"{ROOT}/FreeEnergy/{MOL}/{SYS}/{NGAS}/{KAPPA}/CV_water_gas{INC}/{TARGET}/FE**.gro")[0]
        xtc = glob.glob(f"{ROOT}/FreeEnergy/{MOL}/{SYS}/{NGAS}/{KAPPA}/CV_water_gas{INC}/{TARGET}/FE**.xtc")[0]
    except:
        continue

    u = mda.Universe(gro,xtc)
    zmin = 440 
    zmax = 630

    start_time = i*(N_STEPS)
    ngas = []
    nw = []
    for ts in u.trajectory[-10:]:
        bead = GASES[MOL]
        natoms = u.select_atoms(f"resname {bead} and prop z > {zmin} and prop z < {zmax}").n_atoms
        nwater = u.select_atoms(f"(resname W or resname TW or resname SW) and prop z > {zmin} and prop z < {zmax}").n_atoms
        ngas.append(natoms)
        nw.append(nwater)
    
    mean_ngas = np.mean(ngas)
    std_ngas = np.std(ngas)
    mean_nw = np.mean(nw)
    std_nw = np.std(nw)
    data.append([MOL, mean_ngas, std_ngas, NGAS, TARGET, mean_nw, std_nw])

if data:
    df = pd.DataFrame(data, columns=["mol","av_ngas","std_gas","ngas","coord","av_nwater","std_water"])
    df.to_csv(f"{ROOT}/FreeEnergy/outfiles/ngas_{MOL}_{NGAS}.csv",index=False)
