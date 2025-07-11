
import numpy as np
import matplotlib.pyplot as plt
from .capacity_driven import getCapacityDriven, plot_pile
from .capacity_soils import clay_profile, sand_profile
from .capacity_plots import plot_helical

def getCapacityHelical(profile, soil_type, D, L, d, zlug, Ha, Va, plot=True):
    '''Calculate the vertical and horizontal capacity of a helical pile using a soil profile.
    The calculation is based on the soil profile, anchor geometry and inclined load.

    Parameters
    ----------
    profile : array
        Soil profiles (z, parameters)
            Clay: (z (m), Su (kPa), gamma (kN/m³))
            Sand: (z (m), phi (deg), gamma (kN/m³), Dr (%))
    soil_type : string
        Select soil condition, 'clay' or 'sand'
    D : float 
        Helix diameter (m)
    L : float 
        Shaft length (m)
    d : float 
        Shaft diameter (m)
    zlug : float
        Depth to padeye (m)
    Ha : float
        Horizontal load at pile lug elevation (N)   
    Va : float
        Vertical load at pile lug elevation (N)
    plot : bool
        Plot the p-y curve and the deflection pile condition if True

    Returns
    -------
    Dictionary with capacity, displacements, weight and UC.
    '''

    t = (6.35 + D*20)/1e3            # Helical pile wall thickness (m), API RP2A-WSD
    rhows = 66.90e3                  # Submerged steel specific weight (kN/m3)
    rhow = 10e3                      # Water specific weight (kN/m3) 

    def PileWeight(Len, Dia1, Dia2, tw, rho):
        return ((np.pi/4)*((Dia1**2 - (Dia1 - 2*tw)**2)*Len + (np.pi/4)*Dia2**2*tw))*rho

    if soil_type == 'clay':
        z0, f_Su, f_sigma_v_eff, f_gamma, f_alpha = clay_profile(profile)

        z_helix = np.clip(zlug + (L - D), profile[0, 0], profile[-1, 0])
        Su = f_Su(z_helix)
        sigma_v_eff = max(f_sigma_v_eff(z_helix), 1.0)
        psi_val = Su/sigma_v_eff
        alpha = min(0.5*psi_val**-0.50, 1) if psi_val <= 1.0 else min(0.5 * psi_val**-0.25, 1)

        Nc = min(6.0*(1 + 0.2*d/D), 9)
        Qh = ((np.pi/4)*(D**2 - d**2)*Nc*Su + f_gamma(z_helix)*D)*0.75
        Qs = np.pi*d*L*alpha*Su
        Qu = PileWeight(L, D, d, t, rhows) + Qh + Qs

    elif soil_type == 'sand':
        z0, f_phi, f_sigma_v_eff, f_gamma, f_Dr, f_delta = sand_profile(profile)

        z_helix = np.clip(zlug + (L - D), profile[0, 0], profile[-1, 0])
        gamma = f_gamma(z_helix)
        Dr = f_Dr(z_helix)
        delta = f_delta(z_helix)
        phi = f_phi(z_helix)

        Nq = 0.5*(12*phi)**(phi/54)
        Qh = (np.pi/4)*(D**2 - d**2)*Nq*gamma*z_helix
        Qs = np.pi*d*L*delta*gamma*z_helix
        Qu = PileWeight(L, D, d, t, rhows) + Qh + Qs
        
    # Pile weight (inc. auxilary items) assessed as a factor
    Wp = 1.10*PileWeight(L, D, d, t, (rhows + rhow))    

    # Unity Check based only on vertical capacity
    UC_vertical = Va/Qu

    # Compute horizontal capacity using p-y method
    y, z, results_lateral = getCapacityDriven(profile, soil_type, L, d, zlug, Ha, Va, plot=True)
    if soil_type == 'clay':
        plot_pile(profile_clay, 'clay', y, z, D, L, z0, zlug, hinge_location=None)
    elif soil_type == 'sand':
        plot_pile(profile_sand, 'sand', y, z, D, L, z0, zlug, hinge_location=None)

    Hcap = Ha if 'Plastic moment' not in results_lateral else results_lateral['Bending moment']/abs(zlug) if zlug != 0 else 1e-6
    UC_horizontal = Ha/Hcap if Hcap != 0 else np.inf

    resultsHelical = {
        'Weight pile': Wp,
        'Vertical max.': Qu,
        'Horizontal max.': Hcap,
        'Unity check (vertical)': UC_vertical,
        'Unity check (horizontal)': UC_horizontal,
        'Lateral displacement': results_lateral.get('Lateral displacement', None),
        'Rotational displacement': results_lateral.get('Rotational displacement', None),
        'Bending moment': results_lateral.get('Bending moment', None),
        'Plastic moment': results_lateral.get('Plastic moment', None),
        'Plastic hinge': results_lateral.get('Plastic hinge', None),
        'Hinge location': results_lateral.get('Hinge location', None),
        'p-y model': results_lateral.get('p-y model', None),
        }

    return resultsHelical

if __name__ == '__main__':
    
    profile_clay = np.array([
        [ 1.0, 10, 8.0],
        [ 5.0, 15, 8.5],
        [10.0, 25, 8.5],
        [25.0, 50, 9.0]
    ])

    profile_sand = np.array([
        [ 1.0, 28,  9.5, 40, 60],
        [ 5.0, 34, 10.0, 42, 70],
        [ 8.0, 38, 10.0, 44, 75],
        [15.0, 38, 11.5, 45, 85]
    ])  

    D = 1.8        # Helix diameter (m)
    L = 12.0       # Pile length (m)
    d = 0.8        # Shaft diameter (m)
    zlug = 2       # Padeye depth (m)
    Va = 50e3      # Vertical load (N)
    Ha = 30e3      # Horizontal load (N)

    # === CLAY ===
    # resultsHelical_clay = getCapacityHelical(profile_clay, 'clay', D, L, d, zlug, Ha, Va, plot=True)
    # for key, val in resultsHelical_clay.items():
    #     if isinstance(val, float):
    #         print(f"{key}: {val:.3f}")
    #     else:
    #         print(f"{key}: {val}")
    
    # plot_helical(profile_clay, 'clay', D, L, d, zlug, n_helix=1, spacing=1.0, title='Helical Pile in Clay Profile')

     # === SAND ===
    resultsHelical_sand = getCapacityHelical(profile_sand, 'sand', D, L, d, zlug, Ha, Va, plot=True)
    for key, val in resultsHelical_sand.items():
        if isinstance(val, float):
            print(f"{key}: {val:.3f}")
        else:
            print(f"{key}: {val}")
        
    plot_helical(profile_sand, 'sand', D, L, d, zlug, n_helix=1, spacing=1.0, title='Helical Pile in Sand Profile')
    


