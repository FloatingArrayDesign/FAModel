
import numpy as np
import matplotlib.pyplot as plt
from .support_soils import clay_profile
from .support_plots import plot_plate

def getCapacityPlate(profile_map, location_name, B, L, zlug, beta, Ha, Va, plot=True):
    '''Calculate the plate anchor capacity using clay soil layers from profile_map.
    The calculation is based on the soil profile, anchor geometry and inclined load.

    Parameters
    ----------
    profile_map : list of dict
        Soil profile map with coordinates and layers per location.
    location_name : str
        Name of the location to select the soil profile.
    B : float
        Plate width (m)
    L : float
        Plate length (m)
    zlug : float
        Embedment depth of the main padeye (m)
    beta : float
        Inclination angle of the plate (deg)
    Ha : float
        Applied horizontal load (N)
    Va : float
        Applied vertical load (N)
    plot : bool
        Whether to generate plots.

    Returns
    -------
    Dictionary with Capacity, Weight, UC, etc.
    '''

    # Extract and filter clay layers from profile_map
    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = [layer for layer in profile_entry['layers'] if layer['soil_type'] == 'clay']

    if not layers:
        raise ValueError('Plate anchor capacity model only supports clay soils.')

    # Build the profile array: [[z, Su, gamma], ...]
    profile = []
    for layer in layers:
        profile.append([layer['top'], layer['gamma_top'], layer['Su_top']])
        profile.append([layer['bottom'], layer['gamma_bot'], layer['Su_bot']])
        
    print("layer gamma_top (raw):", layer['gamma_top'])
    print("layer gamma_bot (raw):", layer['gamma_bot'])

    profile = np.array(sorted(profile, key=lambda x: x[0]))

    # Parameters and constants
    Los = 0.05
    B_t = 40
    rhows = 66.90e3  # Submerged steel (N/m3)
    rhow = 10e3      # Seawater (N/m3)

    # Evaluate interpolated Su and gamma
    z0, f_gamma, f_Su, f_sigma_v_eff, f_alpha = clay_profile(profile)
    t = round(B/B_t, 2)
    V_steel = round(B*L*t, 2)
    zlug_B = zlug/B

    # Profile check points
    npts = 10
    z_offsets = np.linspace(-0.5, 0.5, npts)*B*np.sin(np.deg2rad(beta))
    z_points = zlug + z_offsets; print(z_points)

    Su_vals = [f_Su(z) for z in z_points]
    gamma_10 = f_gamma(z_points[2]); print(gamma_10)
    gamma_vals = [f_gamma(z) for z in z_points]; print("gamma_vals:", [f"{val:.2f}" for val in gamma_vals], "N/m3")
    Su = np.mean(Su_vals); print(f"Su: {Su:.2f} Pa")
    gamma = np.mean(gamma_vals); print(f"gamma: {gamma:.2f} N/m3")
    
    print("Profile being sent to clay_profile():")
    for row in profile:
        print(f"z = {row[0]:.2f} m, gamma = {row[1]:.2f} kN/m³, Su = {row[2]:.2f} kPa")

    # Shear strength gradient
    k = np.polyfit(z_points, Su_vals, 1)[0]
    print(f"k: {k:.2f}")

    # Pile weight including auxiliary parts
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

    Nc_final = max(Nco + (gamma*zlug)/Su, Nco_s)
    print(f"Nc_star: {Nco + (gamma*zlug)/Su:.2f}")
    print(f"Nc_star: {Nco_s:.2f}")
    qu = Nc_final*Su
    Tmax = round(qu*(1 - Los)*B*L, 2)
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
    
    return layers, resultsPlate

if __name__ == '__main__':
    profile_map = [
        {
            'name': 'CPT_1',
            'x': 498234, 'y': 5725141,
            'layers': [
                {
                    'top': 0.0, 'bottom': 9.5,
                    'soil_type': 'clay',
                    'gamma_top': 8.0, 'gamma_bot': 8.5,
                    'Su_top': 10, 'Su_bot': 25
                },
                {
                    'top': 9.5, 'bottom': 11.5,
                    'soil_type': 'clay',
                    'gamma_top': 8.5, 'gamma_bot': 8.5,
                    'Su_top': 25, 'Su_bot': 45
                },
                {
                    'top': 11.5, 'bottom': 25.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.5, 'gamma_bot': 9.0,
                    'Su_top': 45, 'Su_bot': 50
                }
            ]
        }
    ]

    B = 2.0
    L = 2.0
    zlug = 10.0
    Ha = 350e3
    Va = 400e3
    alpha = np.rad2deg(np.arctan2(Va, Ha))
    beta = 90 - alpha

    layers, results = getCapacityPlate(profile_map, 'CPT_1', B, L, zlug, beta, Ha, Va)

    print("\n--- Plate Anchor Capacity Results ---")
    for key, val in results.items():
        print(f"{key}: {val:.2f}")
        
    plot_plate(layers, B, L, z0 = layers[0]['top'], zlug=zlug, beta=beta, title='Plate Anchor in Layered Soil')
