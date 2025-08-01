
import numpy as np
from .capacity_driven import getCapacityDriven, plot_pile
from .support_soils import clay_profile, sand_profile
from .support_plots import plot_helical

def getCapacityHelical(profile_map, location_name, D, L, d, zlug, Ha, Va, plot=False):
    '''Calculate the vertical and horizontal capacity of a helical pile using a soil profile.
    The calculation is based on the soil profile, anchor geometry and inclined load.

    Parameters
    ----------
    profile : array
        Soil profiles (z, parameters)
            Clay soil profile (z, Su, gamma)
            Sand soil profile (z, phi, gamma, Dr)
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
        Horizontal load applied at padeye (N)
    Va : float
        Vertical load applied at padeye (N)
    plot : bool
        Plot the p-y curve and the deflection pile condition if True

    Returns
    -------
    y : array
        Lateral displacement at each node (real nodes only)
    z : array
        Node depth positions corresponding to y (m)
    resultsHelical : dict
        Dictionary containing displacements, moment capacity, hinge state and vertical capacity
    '''

    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']

    t = (6.35 + D*20)/1e3            # Helical pile wall thickness (m), API RP2A-WSD
    rhows = 66.90e3                  # Submerged steel specific weight (kN/m3)
    rhow = 10e3                      # Water specific weight (kN/m3) 

    def PileWeight(Len, Dia1, Dia2, tw, rho):
        return ((np.pi/4)*((Dia1**2 - (Dia1 - 2*tw)**2)*Len + (np.pi/4)*Dia2**2*tw))*rho

    z_helix = zlug + (L - D)
    matched_layer = next((layer for layer in layers if layer['top'] <= z_helix <= layer['bottom']), None)
    if matched_layer is None:
        raise ValueError(f"No soil layer found at z = {z_helix:.2f} m")

    if matched_layer['soil_type'] == 'clay':
        profile = [[matched_layer['top'],    matched_layer['Su_top'], matched_layer['gamma_top']],
                   [matched_layer['bottom'], matched_layer['Su_bot'], matched_layer['gamma_bot']]]
        z0, f_Su, f_sigma_v_eff, f_gamma, f_alpha = clay_profile(profile)

        z_helix = np.clip(z_helix, matched_layer['top'], matched_layer['bottom'])
        Su = f_Su(z_helix)
        sigma_v_eff = max(f_sigma_v_eff(z_helix), 1.0)
        psi_val = Su/sigma_v_eff
        alpha = min(0.5*psi_val**-0.50, 1) if psi_val <= 1.0 else min(0.5 * psi_val**-0.25, 1)

        Nc = min(6.0*(1 + 0.2*d/D), 9)
        Qh = ((np.pi/4)*(D**2 - d**2)*Nc*Su + f_gamma(z_helix)*D)*0.75
        Qs = np.pi*d*L*alpha*Su
        Qu = PileWeight(L, D, d, t, rhows) + Qh + Qs

    elif matched_layer['soil_type'] == 'sand':
        profile = [[matched_layer['top'],    matched_layer['phi_top'], matched_layer['gamma_top'], matched_layer['Dr_top']],
                   [matched_layer['bottom'], matched_layer['phi_bot'], matched_layer['gamma_bot'], matched_layer['Dr_bot']]]
        z0, f_phi, f_sigma_v_eff, f_gamma, f_Dr, f_delta = sand_profile(profile)

        z_helix = np.clip(z_helix, matched_layer['top'], matched_layer['bottom'])
        gamma = f_gamma(z_helix)
        Dr = f_Dr(z_helix)
        delta = f_delta(z_helix)
        phi = f_phi(z_helix)

        Nq = 0.5*(12*phi)**(phi/54)
        Qh = (np.pi/4)*(D**2 - d**2)*Nq*gamma*z_helix
        Qs = np.pi*d*L*delta*gamma*z_helix
        Qu = PileWeight(L, D, d, t, rhows) + Qh + Qs
        
    
    Wp = PileWeight(L, D, d, t, (rhows + rhow))    

    # Unity Check based only on vertical capacity
    UC_vertical = Va/Qu

    # Compute horizontal capacity using p-y method
    layers, y, z, results_lateral = getCapacityDriven(profile_map, location_name, D, L, zlug, Ha, Va, plot=False)
    # Plotting
    if plot:
        plot_pile(layers, y, z, D, L, z0=layers[0]['top'], zlug=zlug, hinge_location=None)

    Hcap = results_lateral['Horizontal max.']
    UC_horizontal = Ha/Hcap if Hcap != 0 else np.inf

    resultsHelical = {     
        'Horizontal max.': Hcap,
        'Vertical max.': Qu,
        'Lateral displacement': results_lateral['Lateral displacement'],
        'Rotational displacement': results_lateral['Rotational displacement'],
        'Unity check (horizontal)': UC_horizontal,
        'Unity Check (vertical)': UC_vertical,
        'Weight pile': Wp,}

    if matched_layer['soil_type'] == 'clay':
        resultsHelical['Su @ helix'] = Su
        resultsHelical['alpha'] = alpha
    elif matched_layer['soil_type'] == 'sand':
        resultsHelical['Dr @ helix'] = Dr
        resultsHelical['delta'] = delta
        resultsHelical['phi'] = phi

    return layers, resultsHelical

if __name__ == '__main__':
    
    profile_map = [
        {
            'name': 'CPT_1',
            'x': 498234, 'y': 5725141,
            'layers': [
                {
                    'top': 1.0, 'bottom': 3.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.0, 'gamma_bot': 9.0,
                    'Su_top': 60, 'Su_bot': 50},
                # {
                #     'top': 3.0, 'bottom': 7.0,
                #     'soil_type': 'clay',
                #     'gamma_top': 15.0, 'gamma_bot': 25.0,
                #     'Su_top': 100, 'Su_bot': 150},
                {
                    'top': 3.0, 'bottom': 7.0,
                    'soil_type': 'sand',
                    'gamma_top': 8.0, 'gamma_bot': 8.0,
                    'phi_top': 32, 'phi_bot': 38,
                    'Dr_top': 70, 'Dr_bot': 75},
                {
                    'top': 7.0, 'bottom': 15.0,
                    'soil_type': 'clay',
                    'gamma_top': 25.0, 'gamma_bot': 50.0,
                    'Su_top': 200, 'Su_bot': 400}]
        }
    ]  

    D = 1.5        # Helix diameter (m)
    L = 12.0       # Pile length (m)
    d = 0.5        # Shaft diameter (m)
    zlug = 3       # Padeye depth (m)
    Ha = 30e3      # Horizontal load (N)
    Va = 50e3      # Vertical load (N)

    layers, resultsHelical = getCapacityHelical(profile_map, 'CPT_1', D, L, d, zlug, Ha, Va, plot=True)
    for key, val in resultsHelical.items():
        if isinstance(val, float):
            print(f"{key}: {val:.3f}")
        else:
            print(f"{key}: {val}")
    
    plot_helical(layers, D=D, L=L, d=d, z0=layers[0]['top'], zlug=zlug, n_helix=1, spacing=1.0)
  


