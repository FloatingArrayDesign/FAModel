
import numpy as np
import matplotlib.pyplot as plt
from capacity_soils import clay_profile, sand_profile
from capacity_plots import plot_load

def getTransferLoad(profile, soil_type, Tm, thetam, zlug, line_type, d, w=None, plot=False):
    '''Calculate the transfer load from mudline to main padeye using a layered soil profile.

    Parameters
    ----------
    profile : array
        Soil profile as a 2D array: (z, parameters)
            Clay soil profile (z (m), Su (kPa), gamma (kN/m³))
            Sand soil profile (z (m), phi (deg), gamma (kN/m³), Dr (%))
    soil_type : string
        Select soil condition, 'clay' or 'sand'
    Tm : float
        Mooring line load at mudlevel (N)
    thetam : float
        Mooring line angle at mudlevel (deg)
    zlug : float
        Embedment depth of the lug (m)
    line_type : string
        Select mooring line type, 'chain' or 'wire'
    d : float
        Mooring line diameter (m)
    w : float
        Mooring line unit weight (N/m)
    plot : bool
        Plot the inverse catenary mooring line profile if True

    Returns
    -------
    dict
        Dictionary with transferred load components and depth.
    '''
    
    deltas = 0.2  # discretization step

    # Line mechanical properties
    if line_type == 'chain':
        Et, En = 10, 2.5
    elif line_type == 'wire':
        Et, En = np.pi, 1
    W = w*deltas

    # Soil profile and interpolators
    if soil_type == 'clay':
        Nc = 8.5
        z0, f_Su, _, f_gamma, f_alpha = clay_profile(profile)
    elif soil_type == 'sand':
        nhu = 0.5
        z0, f_phi, _, f_gamma, _, f_delta = sand_profile(profile)

    # Initial values
    T = Tm
    theta = np.deg2rad(thetam)
    drag = 0
    depth = 0.1

    # Tracing lists
    drag_values, depth_values = [], []

    while (zlug - depth) >= 0:
        if soil_type == 'clay':
            Su = f_Su(depth)
            alpha = f_alpha(depth)
            dtheta = (En*d*Nc*Su - W*np.cos(theta))/T*deltas
            dT = (Et*d*alpha*Su + W*np.sin(theta))*deltas
        elif soil_type == 'sand':
            gamma_z = f_gamma(depth)
            delta_z = f_delta(depth)
            phi = f_phi(depth)
            Nq = np.exp(np.pi*np.tan(np.deg2rad(phi)))*(np.tan(np.deg2rad(45 + phi/2)))**2
            dtheta = (En*d*Nq*gamma_z*depth - W*np.cos(theta))/T*deltas
            dT = (Et*d*gamma_z * depth * np.tan(np.deg2rad(delta_z)) + W*np.sin(theta))*deltas

        ddrag  = deltas*np.cos(theta)
        ddepth = deltas*np.sin(theta)

        theta += dtheta
        T -= dT
        drag += ddrag
        depth += ddepth

        if abs(Tm - T) > 0.75*Tm:
            raise Exception(f"Load transfer unrealistic: Tm = {Tm/1e6:.2f} MN vs T = {T/1e6:.2f} MN")

        if not (0 < np.rad2deg(theta) < 90):
            raise Exception(f"Load angle unrealistic: {np.rad2deg(theta):.2f} deg")

        drag_values.append(-drag)
        depth_values.append(-depth)

    Ta = T; thetaa = theta
    # H = Ta*np.cos(thetaa)
    # V = Ta*np.sin(thetaa)

    resultsLoads = {
        'Tm': Tm,
        'thetam': thetam,
        'Ta': Ta,
        'thetaa': thetaa,
        'length': deltas*len(drag_values),
        'drag_values': drag_values,
        'depth_values': depth_values
    }

    return resultsLoads

if __name__ == '__main__':
    
    # Define a clay profile: [depth (m), Su (kPa), gamma (kN/m3)]
    profile_clay = np.array([
        [0.0, 100, 8],
        [2.0, 200, 8.5],
        [5.0, 300, 9],
        [12.0, 400, 9.5]
    ])
    
    # Define a sand profile: [depth (m), phi (deg), gamma (kN/m3), Dr (-)]
    profile_sand = np.array([
        [0.0, 30, 9.5, 60],
        [5.0, 32, 10.0, 70],
        [10.0, 34, 10.2, 80],
        [15.0, 35, 10.5, 90]
    ])
    
    # Input parameters
    soil_type = 'clay'
    Tm = 1.2e7           # Load at mudline (N)
    thetam = 10          # Angle at mudline (deg)
    zlug = 8             # Padeye depth (m)
    line_type = 'chain'
    d = 0.16             # Chain diameter (m)
    w = 4093             # Line weight (N/m)
    
    # Run transfer load calculation
    results = getTransferLoad(profile_clay, soil_type, Tm, thetam, zlug, line_type, d, w, plot=True)
    print("\n--- Transfer Load Results (Clay) ---")
    # for key, val in results.items():
    #     print(f"{key}: {val:.3f}")
        
    plot_load(profile_clay, soil_type, results['drag_values'], results['depth_values'], results['Tm'], results['thetam'], results['Ta'], results['thetaa'], zlug)
      

    # # Input parameters
    # Tm = 1.2e5
    # thetam = 5
    # zlug = 8
    # line_type = 'chain'
    # d = 0.25
    # soil_type = 'sand'
    # w = 4093
    
    # results = getTransferLoad(profile_sand, soil_type, Tm, thetam, zlug, line_type, d, w, plot=True)
    
    # print("\n--- Transfer Load Results (Sand) ---")
    # for key, val in results.items():
    #     print(f"{key}: {val:.3f}")

