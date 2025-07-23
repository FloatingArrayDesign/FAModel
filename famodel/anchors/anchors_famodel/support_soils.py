
import numpy as np
from scipy.interpolate import interp1d

def clay_profile(profile):
    ''' Create interpolated functions for a clay soil profile.
    Calculates Su, effective vertical stress, unit weight and adhesion factor.
    
    Parameters
    ----------
    profile : array
        Clay profile as 2D array: (z, gamma, Su)
        Depth (m), effective unit weight gamma (kN/m³) and undrained shear strength Su (kPa)
    
    Returns
    -------
    z0 : float
        Depth of mudline relative to pile head (m)
    f_gamma : interp1d
        Effective unit weight of the soil, γ'(z) (N/m³)
    f_Su : interp1d
        Undrained shear strength, Su(z) (Pa)
    f_sigma_v_eff : interp1d
        Effective vertical stress, σ'v(z) (Pa)
    f_alpha : function
        Adhesion factor from API correlation, α (-)
    '''

    # Depth of mudline relative to pile head
    z0 = float(profile[0][0])

    # Extract data from profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth = np.concatenate([np.array([z0]),np.array([row[0] for row in profile],dtype=float)]) # (m)
    gamma = np.concatenate([np.array([0]), np.array([row[1] for row in profile],dtype=float)]) # (kN/m3)
    Su    = np.concatenate([np.array([0]), np.array([row[2] for row in profile],dtype=float)]) # (kPa)

    # Calculate sigma_v_eff at each depth
    sigma_v_eff = np.zeros(len(depth))

    for i in range(1, len(depth)):
        sigma_v_eff[i] = sigma_v_eff[i-1] + gamma[i-1]*(depth[i] - depth[i-1])

    # Define interpolation functions
    f_Su          = interp1d(depth, Su*1e3, kind='linear')           # (Pa)
    f_sigma_v_eff = interp1d(depth, sigma_v_eff*1e3, kind='linear')  # (Pa)
    f_gamma       = interp1d(depth, gamma*1e3, kind='linear')        # (N/m3)
    
    # Calculate f_psi and f_alpha at each depth (not as a scalar)
    f_psi = lambda z: f_Su(z)/np.maximum(f_sigma_v_eff(z), 1.0)  
    
    def calc_alpha(psi):
        # Avoid divide-by-zero or log(0) by setting a floor
        psi = np.maximum(psi, 1e-6)
        if np.ndim(psi) == 0:
            psi = float(psi)
            # API-style adhesion factor: two regimes
            return min(0.5*psi**-0.50, 1) if psi <= 1.0 else min(0.5*psi**-0.25, 1)
        else:
            return np.where(
                psi <= 1.0,
                np.minimum(0.5*psi**-0.50, 1),
                np.minimum(0.5*psi**-0.25, 1)
            )

    # Create an interpolated adhesion factor function
    def f_alpha(z):
        psi_val = f_psi(z)
        alpha_val = calc_alpha(psi_val)
        return np.atleast_1d(alpha_val)[0] if np.ndim(alpha_val) == 0 else alpha_val

    return z0, f_gamma, f_Su, f_sigma_v_eff, f_alpha

def sand_profile(profile):
    ''' Create interpolated functions for a sand soil profile.
    Calculates phi, effective stress, unit weight, relative density, and skin friction factor.

    Parameters
    ----------
    profile : array
        Sand profile as 2D array: (z, gamma, phi, Dr)
        Depth (m), effective unit weight, gamma (kN/m³), friction angle, phi (deg) and relative density, Dr (%)

    Returns
    -------
    z0 : float
        Depth of mudline relative to pile head (m)
    f_gamma : interp1d
        Effective unit weight of the soil, γ'(z) (N/m³)
    f_phi : interp1d
        Friction angle, φ(z) (deg)
    f_sigma_v_eff : interp1d
        Effective vertical stress, σ'v(z) (Pa)
    f_Dr : interp1d
        Relative density, Dr(z) (-)
    f_delta : interp1d
        Skin friction factor, δ(z) (-)
    '''

    # Depth of mudline relative to pile head
    z0 = float(profile[0][0])

    # Extract data from profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth = np.concatenate([np.array([z0]),np.array([row[0] for row in profile],dtype=float)]) # (m)
    gamma = np.concatenate([np.array([0]), np.array([row[1] for row in profile],dtype=float)]) # (kN/m3)
    phi   = np.concatenate([np.array([0]), np.array([row[2] for row in profile],dtype=float)]) # (deg)
    Dr    = np.concatenate([np.array([0]), np.array([row[3] for row in profile],dtype=float)]) # (%)
   
    # Calculate sigma_v_eff and static loading factor at each depth
    sigma_v_eff = np.zeros(len(depth))

    for i in range(1, len(depth)):  
        sigma_v_eff[i] = sigma_v_eff[i-1] + gamma[i-1]*(depth[i] - depth[i-1])

    # Define interpolation functions
    f_phi         = interp1d(depth, phi, kind='linear')              # (deg)
    f_sigma_v_eff = interp1d(depth, sigma_v_eff*1e3, kind='linear')  # (Pa)
    f_gamma       = interp1d(depth, gamma*1e3, kind='linear')        # (N/m3)
    f_Dr          = interp1d(depth, Dr, kind='linear')               # (%)
    
    # Define delta as a function of Dr
    def calc_delta(Dr_val):
        if 35 <= Dr_val < 50:
            return 0.29
        elif 50 <= Dr_val < 65:
            return 0.37
        elif 65 <= Dr_val < 85:
            return 0.46
        elif Dr_val >= 85:
            return 0.56
        else:
            return 0  # Default or error value for very low Dr values
        
    # Apply delta calculation to Dr profile
    delta_values = np.array([calc_delta(Dr_val) for Dr_val in Dr])
    f_delta      = interp1d(depth, delta_values, kind='linear')  # Interpolated delta values
    
    return z0, f_gamma, f_phi, f_Dr, f_sigma_v_eff, f_delta
    
def rock_profile(profile):
    ''' Create interpolated functions for a weak rock profile.
    Calculates unconfined compressive strength (UCS) and Young’s modulus (Em).

    Parameters
    ----------
    profile : array
        Rock profile as 2D array: (z, UCS, Em)
        Depth (m), unconfined compressive strength, UCS (MPa), Young's modulus, Em (MPa)

    Returns
    -------
    z0 : float
        Depth of rockline relative to pile head (m)
    f_UCS : interp1d
        Unconfined compressive strength, UCS(z) (Pa)
    f_Em : interp1d
        Young's modulus, Em(z) (Pa)
    '''

    # Depth of rockline relative to pile head
    z0 = float(profile[0][0])

    # Extract data from soil_profile array and zero strength virtual soil layer
    # from the pile head down to the mudline
    depth = np.concatenate([np.array([z0]),np.array([row[0] for row in profile],dtype=float)]) # (m)  
    UCS   = np.concatenate([np.array([0]), np.array([row[1] for row in profile],dtype=float)]) # (MPa)
    Em    = np.concatenate([np.array([0]), np.array([row[2] for row in profile],dtype=float)]) # (MPa)

    # Define interpolation functions
    f_UCS = interp1d(depth, UCS*1e6, kind='linear')  # (Pa)
    f_Em  = interp1d(depth, Em*1e6, kind='linear')   # (Pa)
    
    return z0, f_UCS, f_Em

