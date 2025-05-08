from setup import *
from scipy.signal import argrelextrema

def apply_smoothing(FE):
    FE_smooth = savgol_filter(FE, SMOOTH_WINDOW, SMOOTH_POLY)
    return FE_smooth

def remove_linear_trend(FE, TARGS):
    coeffs = np.polyfit(TARGS[10:-10], FE[10:-10], 1)
    linear_trend = np.polyval(coeffs, TARGS)

    P_lin = coeffs[0] * 4.11433402
    FE_corrected = FE - linear_trend
    FE_corrected = FE_corrected - FE_corrected.min() 

    return FE_corrected, P_lin

def calc_intrusion_pressure_from_barrier(FE, TARGS, minima_idx, maxima_idx):
    # Last minimum and the last maximum
    try:
        last_min_idx = minima_idx[0]
        last_max_idx = maxima_idx[0]
    except:
        return None

    FE_min_last = FE[last_min_idx]
    FE_max_last = FE[last_max_idx]

    deltaG_barrier = FE_max_last - FE_min_last  # in kBT

    # Estimate delta_V
    deltaV = TARGS[last_max_idx] - TARGS[last_min_idx] # in nm^3
    Pint = deltaG_barrier / deltaV * 4.14  # in MPa
    return Pint

def calc_extrusion_pressure_from_barrier(FE, TARGS, minima_idx, maxima_idx):
    # Last minimum and the last maximum
    try:
        last_min_idx = minima_idx[-1]
        last_max_idx = maxima_idx[-1]
    except:
        return None

    FE_min_last = FE[last_min_idx]
    FE_max_last = FE[last_max_idx]

    deltaG_barrier = - (FE_max_last - FE_min_last)  # in kBT

    # Estimate delta_V
    deltaV = TARGS[last_min_idx] - TARGS[last_max_idx] # in nm^3
    Pext = deltaG_barrier / deltaV * 4.14  # in MPa
    return Pext


def calc_inflection_points(FE, TARGS):
    dFE = np.gradient(FE, TARGS)  # First derivative
    dFE_smooth = savgol_filter(dFE, SMOOTH_WINDOW, SMOOTH_POLY)  # Smooth dFE

    # Find Inflection Points (zero-crossings of d²FE/dx²)
    d2_FE = np.gradient(dFE_smooth, TARGS)  # Second derivative
    inflection_indices = np.where(np.diff(np.sign(d2_FE)))[0]

    if len(inflection_indices) < 2:
        print(f"Skipping {MOL} ({GAS}, NGAS={NGAS}) - Less than 2 inflection points found.")
        return None, None, None
    elif len(inflection_indices) > 2:
        int_idx = [inflection_indices[0],inflection_indices[2]]
    else:
        int_idx = [inflection_indices[0]]

    # Extract first and last inflection points
    ext_idx = inflection_indices[-1]

    # Use np.concatenate instead of unpacking into a list
    inflection_coords = np.concatenate((TARGS[int_idx], [TARGS[ext_idx]]))
    inflection_points = np.concatenate((dFE_smooth[int_idx], [dFE_smooth[ext_idx]]))
    inflection_FE = np.concatenate((FE_corrected_smooth[int_idx], [FE_corrected_smooth[ext_idx]]))
    return inflection_coords, inflection_points, inflection_FE, 

def plot_inflection_points(TARGS, FE_corrected, FE_corrected_smooth, dFE, dFE_smooth, inflection_coords, inflection_FE, inflection_points):
        plt.figure(figsize=(5, 3))
        plt.plot(TARGS, FE_corrected, label="FE_corrected (Raw)", color='gray', linestyle="dashed", alpha=0.6)
        plt.plot(TARGS, FE_corrected_smooth, label="FE_corrected (Smoothed)", color='b')

        for i, coord in enumerate(inflection_coords):
            plt.scatter(coord, inflection_FE[i], c='r', label=f"Inflection Point {i+1}", zorder=3)
        plt.xlabel("Coord")
        plt.ylabel("FE_corrected (kbT)")
        plt.legend()
        plt.title(f"Free Energy Profile - {GAS}, NGAS={NGAS}")
        plt.show()

        plt.figure(figsize=(8, 5))
        plt.plot(TARGS, dFE, label="dFE/dx (Raw)", color='gray', linestyle="dashed", alpha=0.6)
        plt.plot(TARGS, dFE_smooth, label="dFE/dx (Smoothed)", color='g')
        for i, coord in enumerate(inflection_coords):
            plt.scatter(coord, inflection_points[i], c='r', label=f"Inflection Slope {i+1}", zorder=3)
        plt.xlabel("Coord")
        plt.ylabel("dFE/dx")
        plt.legend()
        plt.title(f"First Derivative - {GAS}, NGAS={NGAS}")
        plt.show()

def calc_intrusion_pressure_from_inflection(inflection_points, inflection_coords, inflection_FE):
    ipoints_pint = inflection_points[0:2].copy()
    icoords_pint = inflection_coords[0:2].copy()
    iFE_pint = inflection_FE[0:2].copy()

    PINT = np.array(inflection_points * 4.11433402)[0]
    return PINT, icoords_pint, ipoints_pint, iFE_pint

def calc_extrusion_pressure_from_inflection(inflection_points, inflection_coords, inflection_FE):
    ipoints_pext = inflection_points[-1].item()
    icoords_pext = inflection_coords[-1].item()
    iFE_pext = inflection_FE[-1].item()

    PEXT = np.array(inflection_points * 4.11433402)[-1]
    return PEXT, icoords_pext, ipoints_pext, iFE_pext

if __name__ == "__main__":

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.integrate import cumulative_trapezoid as cumtrapz
    from scipy.signal import savgol_filter  # Import Savitzky-Golay filter

    # ----- Load Data -----
    ROOT = "/Users/frasera/Ricerca"
    SYS = "C8,1.2"
    INC = "_inc_50"
    KAPPA = 0.01

    SHIFT_RANGE = 1  # Number of indices to shift left and right
    SMOOTH_WINDOW = 7  # Window size for smoothing (must be odd)
    SMOOTH_POLY = 2  # Polynomial order for Savitzky-Golay filter

    int_inf_points, ext_inf_points = [], []
    Pint_Pext_list = []
    dfs = []

    print(MOLS)
    for MOL in MOLS:
        GAS = GASES[MOL]

        for NGAS in [0,10,100, 1000,3000,5000]:
            try:
                df = pd.read_csv(f"{ROOT}/FreeEnergy/outfiles/meanF_{KAPPA}_{SYS}_{GAS}_{NGAS}.csv")
            except Exception as err:
                continue

            # ----- Calculate Free Energy (FE) -----
            FE = cumtrapz(y=df["meanF"].values, x=df["coord"].values) / 2.430527019  # in kbT
            FE = np.insert(FE, 0, 0)  # Match length with coord
            TARGS = df["coord"].values * float(2.99e-2)

            # ----- Remove Linear Trend -----
            FE_corrected, P_lin = remove_linear_trend(FE, TARGS)

            # ----- Apply Smoothing -----
            FE_corrected_smooth = apply_smoothing(FE_corrected)

            # ---- Find minima and maxima -----
            minima_idx = argrelextrema(FE_corrected_smooth, np.less)[0]
            maxima_idx = argrelextrema(FE_corrected_smooth, np.greater)[0]

            # ----- Calculate Derivatives -----
            FE_min = FE_corrected_smooth[minima_idx[0]]
            FE_corrected_smooth -= FE_min
            inflection_coords, inflection_points, inflection_FE,  = calc_inflection_points(FE_corrected_smooth, TARGS)
        
            # ----- Calculate Intrusion and Extrusion Pressures -----
            PINT, icoords_pint, ipoints_pint, iFE_pint = calc_intrusion_pressure_from_inflection(inflection_points, inflection_coords, inflection_FE)
            PEXT, icoords_pext, ipoints_pext, iFE_pext = calc_extrusion_pressure_from_inflection(inflection_points, inflection_coords, inflection_FE)
            PINT_barrier = calc_intrusion_pressure_from_barrier(FE_corrected_smooth, TARGS, minima_idx, maxima_idx)
            PEXT_barrier = calc_extrusion_pressure_from_barrier(FE_corrected_smooth, TARGS, minima_idx, maxima_idx)

            #plot_points(TARGS, FE_corrected, FE_corrected_smooth, dFE, dFE_smooth, inflection_coords, inflection_FE, inflection_points)

            print("#-----------------")
            print(f"{MOL} ({GAS}, NGAS={NGAS}):")
            print(f"Pint (inflection): {PINT} / Pext (inflection): {PEXT}")
            print(f"Pint (barrier): {PINT_barrier} / Pext (barrier): {PEXT_barrier}")

            try:
                Delta_FE = FE_corrected_smooth[maxima_idx[-1]] - FE_corrected_smooth[minima_idx[-1]]
            except:
                Delta_FE = None

            Pint_Pext_list.append([MOL.replace("_"," "), NGAS, PINT, PEXT, PINT_barrier, PEXT_barrier, NGAS/(2*20**3), Delta_FE])

            df = pd.DataFrame({'coord':TARGS,
                            'FE':FE_corrected_smooth, 
                            'mol':[MOL.replace("_"," ")]*len(TARGS),
                            'ngas':[NGAS]*len(TARGS),
                            'C':[NGAS/(2*20**3)]*len(TARGS)
                            })
            dfs.append(df)

            # ----- Create DataFrames for Inflection Points -----
            vals_int = np.array([*icoords_pint, *iFE_pint,
                                MOL.replace("_"," "), 
                                NGAS, NGAS/(2*20**3)])
            df_int_inf_pts = pd.DataFrame([vals_int],
                                          columns=['x_int','x_int_1',
                                                   'y_int','y_int_1',
                                                    'mol','ngas','C'])
            int_inf_points.append(df_int_inf_pts)

            
            vals_ext = np.array([icoords_pext, iFE_pext,
                                MOL.replace("_"," "), 
                                NGAS, NGAS/(2*20**3)])
            df_ext_inf_pts = pd.DataFrame([vals_ext], 
                                          columns=['x_ext','y_ext',
                                                    'mol','ngas','C'])
            ext_inf_points.append(df_ext_inf_pts)

    # Save all dataframes to CSV
    df_Pint_Pext = pd.DataFrame(Pint_Pext_list, columns=['mol','ngas','Pint_inflection','Pext_inflection',
                                                         'Pint_barrier','Pext_barrier','C','Delta_FE'])
    df_Pint_Pext.to_csv(f"{ROOT}/FreeEnergy/outfiles/Pint_Pext_all_{KAPPA}_{SYS}.csv", index=False)

    FEs = pd.concat(dfs)
    FEs.to_csv(f'{ROOT}/FreeEnergy/outfiles/FE_all.csv',index=False)

    InfPts= pd.merge(pd.concat(int_inf_points), pd.concat(ext_inf_points), on=['mol','ngas','C'])
    InfPts.to_csv(f'{ROOT}/FreeEnergy/outfiles/inf_points_all.csv',index=False)