import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import MDAnalysis as mda
from tqdm.notebook import tqdm
import glob
from natsort import natsorted
from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

ROOT = "/scratch/rasera"

MOLS = ["Chloroform","N2","TN6","Diethyl_ether","Oxygen","Ethylene","No_gas","Fullerene","Bare","SN2","SX1","SC4","TC1"]

GASES = {"Chloroform":"CLF",
         "Diethyl_ether":"ETH",
         "Ethylene":"ETY",
         "No_gas":"No_gas",
         "TN6":"NI2",
         "Fullerene":"F15",
         "N2":"XEN",
         "Oxygen":"OXY",
         "Bare":"Bare",
         "SN2":"SN2",
         "SX1":"SX1",
         "SC4":"SC4",
         "TC1":"TC1"}

BEADS = {"Chloroform":"CX",
         "Diethyl_ether":"CO",
         "Ethylene":"C4",
         "No_gas":"No_gas",
         "TN6":"TN6",
         "Fullerene":"F15",
         "Bare":"",
         "N2":"N2",
         "Oxygen":"TC3",
         "SN2":"SN2",
         "SX1":"SX1",
         "SC4":"SC4",
         "TC1":"TC1"}

CONC = {"No gas":0,
        "Chloroform":0.040356,
        "Diethyl ether":0.558965,
        "Ethylene":0.420741,
        "Nitrogen":0.000408,
        "Xenon":0.0275}

EXP_PINT = {
    20: {
        "Nitrogen": 17.9,
        "Xenon": 17.1,
        "Diethyl ether": 17.8,
        "Ethylene": 17.9,
        "Ar": 18.9
    },
    25: {
        "Nitrogen": 17.84,
        "Xenon": 17.7,
        "Diethyl ether": 17.8,
        "Ethylene": None,  # Not measured at 25°C
        "Ar": None         # Not measured at 25°C
    },
    50: {
        "Nitrogen": 16.6,
        "Xenon": 17.5,
        "Diethyl ether": 17.4,
        "Ethylene": None,  # Not measured at 50°C
        "Ar": None         # Not measured at 50°C
    },
    80: {
        "Nitrogen": 16.5,
        "Xenon": 15.6,
        "Diethyl ether": 16.2,
        "Ethylene": None,  # Not measured at 80°C
        "Ar": None         # Not measured at 80°C
    }
}

EXP_PEXT = {
    20: {
        "Nitrogen": 1.38,
        "Xenon": 1.98,
        "Diethyl ether": 1.84,
        "Ethylene": 1.52,
        "Ar": 1.1
    },
    25: {
        "Nitrogen": 1.3,
        "Xenon": 1.9,
        "Diethyl ether": 2.1,
        "Ethylene": None,  # Not measured at 25°C
        "Ar": None         # Not measured at 25°C
    },
    50: {
        "Nitrogen": 1.67,
        "Xenon": 2.2,
        "Diethyl ether": 2.2,
        "Ethylene": None,  # Not measured at 50°C
        "Ar": None         # Not measured at 50°C
    },
    80: {
        "Nitrogen": 2.4,
        "Xenon": 3.1,
        "Diethyl ether": 3.1,
        "Ethylene": None,  # Not measured at 80°C
        "Ar": None         # Not measured at 80°C
    }
}

MOL_MAP = {
    'TN6': 'TN6', 
    'Chloroform': 'CHCl$_3$', 
    'Ethylene': 'C$_2$H$_4$', 
    'Diethyl ether': '(C$_2$H$_5$)$_2$O',
    'N2': 'N2',
    'Oxygen': 'O$_2$',
    'SN2': 'SN2',
    'SX1': 'SX1',
    'SC4': 'SC4',
    'TC1': 'TC1'
}

palette = ["#EF476F", "#FFD166", "#06D6A0", "#118AB2", "#073B4C"]
