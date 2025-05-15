
import numpy as np
import matplotlib.pyplot as plt
from .capacity_soils import clay_profile
from .capacity_plots import plot_plate

def getCapacityPlate(profile, soil_type, B, L, zlug, beta, Ha, Va, plot=True):
    '''Calculate the plate anchor capacity using a full clay soil profile and return capacity + UC.
    The calculation is based on the soil profile, anchor geometry and inclined load.

    Parameters
    ----------
    profile : ndarray
        Clay soil profile as a 2D array: [[depth, Su, gamma], ...]
    soil_type : str
        Currently only 'clay' is supported.
    B : float
        Plate width (m)
    L : float
        Plate length (m)
    zlug : float
        Embedment depth of the main padeye (m)
    beta : float
        Inclination angle of the plate (deg)
    Ha : float
        Applied horizontal load (kN)
    Va : float
        Applied vertical load (kN)
    plot : bool
        Placeholder for future use.

    Returns
    -------
    Dictionary with capacity, weight and UC.
    '''
    
    Los = 0.05
    B_t = 40
    rhows = 66.90e3        # Submerged steel specific weight (kN/m3)
    rhow = 10e3            # Water specific weight (kN/m3) 
    
    # Extract soil parameters from profile
    z0, f_Su, _, f_gamma, _ = clay_profile(profile)

    # Geometry
    t = round(B/B_t, 2)
    V_steel = round(B*L*t, 2)
    zlug_B = zlug/B
    
    # Define 5 evaluation points along inclined width
    N = 5
    z_offsets = np.linspace(-0.5, 0.5, N)*B*np.sin(np.deg2rad(beta))
    z_points = zlug + z_offsets; print(z_points)

    # Evaluate Su and gamma at these points
    Su_vals = [f_Su(z) for z in z_points]
    gamma_10 = f_gamma(z_points[2]); print(gamma_10)
    gamma_vals = [f_gamma(z) for z in z_points]; print("gamma_vals:", [f"{val:.2f}" for val in gamma_vals], "N/m3")
    Su = np.mean(Su_vals); print(f"Su: {Su:.2f} Pa")
    gamma = np.mean(gamma_vals); print(f"gamma: {gamma:.2f} N/m3")
    
    print("Profile being sent to clay_profile():")
    for row in profile:
        print(f"z = {row[0]:.2f} m, Su = {row[1]:.2f} kPa, gamma = {row[2]:.2f} kN/m³")

    # Compute shear strength gradient k from linear fit
    k = np.polyfit(z_points, Su_vals, 1)[0]
    print(f"k: {k:.2f}")
    
    # Pile weight (inc. auxiliary elements) assessed as a factor
    Wp = 1.35*V_steel*(rhows + rhow)

    # Capacity factors
    Nco_0_0  = 2.483*np.log(zlug_B) + 1.974
    Nco_90_0 = 2.174*np.log(zlug_B) + 3.391
    kBSh = k*B/Su
    print(f"kBSh: {kBSh:.2f}")

    f0  = np.where(zlug_B < 4, 1.77*(zlug_B**0.3) - 1.289, 0.192*zlug_B + 0.644)
    f90 = np.where(zlug_B < 4, 0.68*(zlug_B**0.5) - 0.410, 0.153*zlug_B + 0.341)

    S_kB_0  = 1 - f0*kBSh
    S_kB_90 = 1 - f90*kBSh
    Nco_0  = S_kB_0*Nco_0_0
    Nco_90 = S_kB_90*Nco_90_0
    Nco = Nco_0 + (Nco_90 - Nco_0)*(beta/90)**2

    Nco_s_0_0  = np.where(2.90*zlug_B + 6.02 <= 11.59, 2.90*zlug_B + 6.02, 11.596)
    Nco_s_90_0 = np.where(2.72*zlug_B + 4.02 <= 11.59, 2.72*zlug_B + 4.02, 11.596)

    S_s_kB_0 = np.where(zlug_B <= 2, 1 + (0.8 - 0.3*zlug_B)*kBSh - (0.383*kBSh**1.36), 1)
    f90s = np.where(zlug_B <= 3, 0.267*zlug_B, 0.6)
    S_s_kB_90 = 1 - f90s*kBSh
    Nco_s_0 = S_s_kB_0*Nco_s_0_0
    Nco_s_90 = S_s_kB_90*Nco_s_90_0
    Nco_s = Nco_s_90 + (Nco_s_0 - Nco_s_90)*((90 - beta)/90)**2

    # Existing capacity factor and base pressure
    Nc_final = max(Nco + (gamma*zlug)/Su, Nco_s) # anchor pullout capacity factor [kN]
    print(f"Nc_star: {Nco + (gamma*zlug)/Su:.2f}")
    print(f"Nc_star: {Nco_s:.2f}")
    qu = Nc_final*Su                   # Bearing pressure capacity of the anchor plate
    Tmax = round(qu*(1 - Los)*B*L, 2)  # Bearing tension force capacity of the anchor plate
    Hmax = Tmax*np.cos(np.deg2rad(90 - beta))
    Vmax = Tmax*np.sin(np.deg2rad(90 - beta))
    
    Ta = np.sqrt(Ha**2 + Va**2)
    UC = Ta/Tmax

    resultsPlate = {
        'Capacity': Tmax,
        'Horizontal max.': Hmax,
        'Vertical max.': Vmax,
        'Unity check': UC,
        'Weight plate': Wp
    }
    
    return resultsPlate

if __name__ == '__main__':
    
    profile_clay = np.array([
        [ 0.0, 10, 8.0],
        [ 9.5, 25, 8.5],
        [11.5, 45, 8.5],
        [25.0, 50, 9.0]
    ])
    
    B = 2.0                                 # Plate width (m)
    L = 2.0                                 # Plate length (m)
    zlug = 10.0                             # Padeye depth (m)
    Ha = 350e3                              # Horizontal load (N)
    Va = 400e3                              # Vertical load (N)
    alpha = np.rad2deg(np.arctan2(Va, Ha))  # Load angle from horizontal (deg)
    beta = 90 - alpha                       # Plate angle after keying (m)
    
    results = getCapacityPlate(profile_clay, 'clay', B, L, zlug, beta, Ha, Va)
    print("\n--- Plate Anchor Capacity Results ---")
    for key, val in results.items():
        print(f"{key}: {val:.2f}")
        
    # plot_plate(profile_clay, 'clay', B, L, zlug, beta, title='Inclined Plate Anchor in Clay')
    
