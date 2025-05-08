from setup import *

def extrapolate_FE(df, n_fit=2, n_extrap=4):  # Separate parameters for fit and extrapolation
    extrapolated_rows = []

    for mol in df['mol'].unique():
        df_mol = df[df['mol'] == mol]
        for ngas_val in df_mol['ngas'].unique():
            subdf = df_mol[df_mol['ngas'] == ngas_val].sort_values('coord')

            if len(subdf) < n_fit + n_extrap:
                continue

            # Fit line to first and last n_fit points
            start = subdf.iloc[:n_fit]
            end = subdf.iloc[-n_fit:]

            slope_start, intercept_start = np.polyfit(start['coord'], start['FE'], 1)
            slope_end, intercept_end = np.polyfit(end['coord'], end['FE'], 1)

            delta_start = start['coord'].iloc[1] - start['coord'].iloc[0]
            delta_end = end['coord'].iloc[-1] - end['coord'].iloc[-2]

            # Extrapolate n_extrap points before and after the data range
            coords_start = [start['coord'].min() - delta_start * (i + 1) for i in reversed(range(int((255-start['coord'].min())//delta_start)))]
            coords_end = [end['coord'].max() + delta_end * (i + 1) for i in range(int((325-end['coord'].max())//delta_end))]

            FE_start = [slope_start * c + intercept_start for c in coords_start]
            FE_end = [slope_end * c + intercept_end for c in coords_end]

            C_val = subdf['C'].iloc[0] if 'C' in subdf.columns else None

            # Collect rows for extrapolated points
            extrapolated_rows.extend([
                {'mol': mol, 'ngas': ngas_val, 'coord': c, 'FE': f, 'C': C_val}
                for c, f in zip(coords_start, FE_start)
            ])
            extrapolated_rows.extend([
                {'mol': mol, 'ngas': ngas_val, 'coord': c, 'FE': f, 'C': C_val}
                for c, f in zip(coords_end, FE_end)
            ])

    # Create the DataFrame with extrapolated rows and append to original dataframe
    df_extra = pd.DataFrame(extrapolated_rows)
    df = pd.concat([df, df_extra], ignore_index=True)
    df = df.sort_values(['mol', 'ngas', 'coord'])

    return df

def coord_num(rij,r0):
    d0 = 0.0
    n = 6
    m = 12

    sij = ( 1 - ((rij - d0)/r0)**n )/( 1 - ((rij - d0)/r0)**m )
    return sij

from PIL import Image
def make_gif(frame_folder,out_name):
    frames = [Image.open(image) for image in natsorted(glob.glob(f"{frame_folder}/*.PNG"))]
    frame_one = frames[0]
    frame_one.save(f"{out_name}", format="GIF", append_images=frames,
               save_all=True, duration=10, loop=10)
    
def plot_FE(mean_force, std_force, coord, mol, color="k", ax=None, cut=None, transform_coord=True, label="default", ls="-"):
    MEANF, COORD, STDF = mean_force, coord, std_force
    FE = cumtrapz(y=MEANF, x=COORD) / 2.430527019   # Convert to k_BT
    FE = np.insert(FE, 0, 0)  # Ensure same length as COORD
    
    dx = np.diff(COORD)  
    dx = np.insert(dx, 0, dx[0])
    FE_STD = np.sqrt(np.cumsum((STDF * dx / 2.430527019) ** 2))

    FE_STD = np.insert(FE_STD, 0, 0) 
    if transform_coord:
        xlabel = "Volume, $nm^3$"
    else:
        xlabel = "Coord"

    TARGS = COORD * 2.99e-2

    # --- Remove Linear Slope ---
    coeffs = np.polyfit(TARGS, FE, 1)  # Linear fit
    linear_trend = np.polyval(coeffs, TARGS)
    FE_corrected = FE - linear_trend  # Subtract linear trend

    # --- Find First Minimum and First Maximum ---
    local_minima = argrelextrema(FE_corrected, np.less)[0]  # Indices of local minima
    local_maxima = argrelextrema(FE_corrected, np.greater)[0]  # Indices of local maxima

    if len(local_minima) == 0 or len(local_maxima) == 0:
        print(f"No local minima or maxima found for {mol}. Skipping...")
        PINT = 0
        PEXT = np.inf
    else:
        # Find the first minimum
        first_min_idx = local_minima[0]
        first_min_FE = FE_corrected[first_min_idx-1]
        first_min_TARG = TARGS[first_min_idx-1]

        second_min_idx = local_minima[1]
        second_min_FE = FE_corrected[second_min_idx-1]
        second_min_TARG = TARGS[second_min_idx-1]

        # Find the first maximum after the first minimum
        first_max_idx = local_maxima[local_maxima > first_min_idx][0]  # First maximum after the first minimum
        first_max_FE = FE_corrected[first_max_idx-1]
        first_max_TARG = TARGS[first_max_idx-1]

        second_max_idx = local_maxima[(local_maxima < second_min_idx)][-1]  # First maximum after the first maximum and before second minimum
        second_max_FE = FE_corrected[second_max_idx-1]
        second_max_TARG = TARGS[second_max_idx-1]

        # --- Compute Slope Between First Minimum and First Maximum ---
        delta_FE = first_max_FE - first_min_FE
        delta_TARG = first_max_TARG - first_min_TARG
        force_slope = delta_FE / delta_TARG  # Slope in k_B T / nm^3

        delta_FE_ext = second_max_FE - second_min_FE
        delta_TARG_ext = second_max_TARG - second_min_TARG
        force_slope_ext = delta_FE_ext / delta_TARG_ext  # Slope in k_B T / nm^3

        # --- Estimate Pressure ---
        PINT = force_slope * 4.11433402 # Convert slope to pressure in MPa
        PEXT = force_slope_ext * 4.11433402

        print(f"{mol}: First minimum at {first_min_TARG:.3f} nm^3, First maximum at {first_max_TARG:.3f} nm^3")
        print(f"Estimated P_int = {PINT:.2f} MPa")
        print(f"{mol}: Second minimum at {second_min_TARG:.3f} nm^3, First maximum at {second_max_TARG:.3f} nm^3")
        print(f"Estimated P_int = {PEXT:.2f} MPa")

    # --- Smooth Interpolation ---
    if transform_coord:
        X = TARGS
        X = X - np.min(X)
    else:
        X = COORD

    Y = FE_corrected
    YSTD = STDF
    
    if cut:
        Y = Y[slice(cut[0],cut[1])]
        X = X[slice(cut[0],cut[1])]
        YSTD = YSTD[slice(cut[0],cut[1])]

    #Y = Y - Y[0]
    # --- Plot --- #
    if ax:
        x_smooth = np.linspace(np.min(X), np.max(X), 100)

        cubic_interp_FE = interp1d(X, Y, kind='cubic')
        y_smooth_cubic = cubic_interp_FE(x_smooth)

        cubic_interp_STD = interp1d(X, YSTD, kind='cubic')
        y_smooth_std = cubic_interp_STD(x_smooth)

        MOLNAME = mol.replace("_", " ")
        if label == "default":
            label = f"{MOLNAME}, $P_{{int}}={PINT:.2f}$ MPa"
        elif not label:
            label = ""
            
        ax.plot(x_smooth, y_smooth_cubic, c=color, lw=3, zorder=-1, label=label, ls=ls)
        #ax.fill_between(x_smooth, y_smooth_cubic - y_smooth_std, y_smooth_cubic + y_smooth_std, 
        #                color=color, alpha=1, zorder=-2)

        #ax.scatter(X, Y, fc="white", ec=color, s=15, zorder=1)
        if label:
            ax.legend(frameon=False)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("FE, $k_BT$")

    return X, Y, YSTD, PINT, PEXT

def quick_plot(file, gas=None):
    u = mda.Universe(file)

    pis1 = u.select_atoms("resname PIS").positions
    pis2 = u.select_atoms("resname PIS1").positions
    bulk = u.select_atoms("resname BULK").positions
    water = u.select_atoms("resname W or resname TW or resname SW").positions

    plt.scatter(pis1[:,0],pis1[:,2],s=1,c="k")
    plt.scatter(pis2[:,0],pis2[:,2],s=1,c="k")
    plt.scatter(bulk[:,0],bulk[:,2],s=1,c="k")
    plt.scatter(water[:,0],water[:,2],s=0.1,c="b")

    minz, maxz = 440, 630
    plt.hlines(minz,-10,210)
    plt.hlines(maxz,-10,210)


    if gas is not None:
        gas_mols = u.select_atoms(f"resname {gas}")
        gas_pos = gas_mols.positions
        plt.scatter(gas_pos[:,0],gas_pos[:,2],s=5,c="r")

        gas_inside = gas_mols.select_atoms(f'prop z > {minz} and prop z < {maxz}')
        print(f'ngas = {gas_inside.n_atoms}')

if __name__ == "__main__":
    import sys as system

    ROOT = "/scratch/rasera"
    MOL = system.argv[1]
    GAS = system.argv[2]
    NGAS = system.argv[3]
    SYS = "C8,1.2"
    INC = "_inc_50"
    KAPPA = 0.01

    REST_TIME = 50000; EQ_TIME = 500000
    N_STEPS = REST_TIME + EQ_TIME; END_STEPS = 550000
    AT0 = 8400; ATF = 11000; INCREMENT_AT = 50

    TARGETS = np.arange(AT0, ATF+INCREMENT_AT, INCREMENT_AT)

    fig, ax = plt.subplots(3,1, figsize=(10,10), sharex=False, gridspec_kw={"hspace":0.3})

    ax[0].set_ylabel("Force, kJ/mol")
    ax[0].set_xlabel("time, ns")
    ax[1].set_ylabel("Prob. dens.")
    ax[1].set_xlabel("$n_{coord}$")

    ax[0].axvline(x=140, color='r', linestyle='--')
    ax[0].axvline(x=150, color='r', linestyle='--')
    ax[0].text(145, -2, "10 ns", ha="center", bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.3'))
    ax[0].text(0.9, 0.1, f"$\kappa=${KAPPA}", ha="right", transform=ax[0].transAxes)

    OUT = []
    STEP = 0

    for i,TARGET in enumerate(tqdm(TARGETS)):
        INIT_TIME = STEP*0.02
        STEP += EQ_TIME

        try:
            CV = pd.read_csv(f"{ROOT}/FreeEnergy/Simulations/{MOL}/{SYS}/{NGAS}/{KAPPA}/CV_water_gas{INC}/{TARGET}/cv.dat", on_bad_lines="skip", 
                        header=None, skiprows=1, usecols=[0,1,2], names=["time","coord","bias"], sep="\s+", dtype={"time":float,"coord":float,"bias":float})
        except Exception as e:
            continue
        
        if len(CV) < 10000:
            continue
        
        cut = CV[(CV["time"]>(REST_TIME*0.02))]
        diff = cut["coord"].values - TARGETS[i]
        time = (cut["time"].values + INIT_TIME)/1000
        coord = cut["coord"].values
        force = -KAPPA*diff

        idx = np.where(np.abs(force) < 10)
        ax[0].plot(time[idx], force[idx])
        ax[1].hist(coord[idx], density=True, histtype="stepfilled", alpha=0.6)

        OUT.append([time[idx].mean(), TARGETS[i], np.mean(force[idx]), np.std(force[idx])])

    dfout = pd.DataFrame(np.array(OUT), columns=["time","coord","meanF","stdF"])
    dfout.to_csv(f"{ROOT}/FreeEnergy/outfiles/meanF_{KAPPA}_{SYS}_{GAS}_{NGAS}.csv", index=False)

    X, Y, YSTD, PINT, PEXT = plot_FE(dfout["meanF"], dfout["stdF"], dfout["coord"], MOL, "xkcd:tomato", ax=ax[2])

    fig.savefig(f"{ROOT}/FreeEnergy/plots/Forces_{KAPPA}_{SYS}_{GAS}_{NGAS}.png", dpi=250)
