"""Plate anchor capacity calculation functions in 'level 2' model, designed
especially for SEPLA (Suction Embedded Plate Anchor), currently set up for
computing the load perpendicular to the anchor plate.
Lead author: Ahmed Radwan. 
"""

import numpy as np

def getCapacityPlate(A, B_t_aspect=40, Hs=20, Bita=30, Los=0.05, 
                     soil_type='clay', gamma=0, So0=2.39, k=1.41):
    '''Calculate the inclined load capacity of a Suction Caisson Anchor in sand or clay.
    The calculation is based on the soil properties, anchor geometry, and the angle of inclined load.  
    
    Parameters
    ----------
    A : float 
        Plate area, assumed to be square so that width B = sqrt(A). [m^2]
    B_t_aspect : float 
        The aspect ratio of the anchor width and its thickness, default is 40.
    soil_type : string
        Specify 'sand' or 'clay'. This affects what other soil parameters are used.
    gamma: float 
        The effective unit weight of the sand soil, default is 9 for sand, 4.7 for clay. [kN/m3] 
    Su0 : float 
        The Undrained shear strength at the mudline (clay only), default is 2.39. [kPa]
    k : float 
        The Undrained shear strength gradient (clay only), default is 1.41. [kPa/m]
    Bita : float
        The Angle of the Plate after Keying Process, default 30. [deg]
    Hs   : float
        The Specified of the Plate embedded depth, default 20. [m]
    Los  : float
        The undrained shear strength lost fraction due to the keying process, default 0.05.
    
    Returns
    -------
    
    '''

    # >>> is gamma with or without water?? <<<

 
    St_density= 7850.0 # The steel density in kg/m3
 
    B = np.sqrt(A)    # Anchor width (and length, approximated as square for now) [m]
    H_B = Hs/B        # Anchor depth range/ The width of the Plate 
    t = B/B_t_aspect    # thickness of the plate, which it depends on the width [m]
    
    #t=np.sqrt(0.006*A)/4  
    V = A*t           # steel volume [m3]
    W = V*St_density  # plate weight [kg]

    Suh = Su0 + k * Hs     # The final clay undrained shear strength


    # ----- anchor pullout capacity -----

    # Anchor Pullout capacity factor in weightless clay with breakaway base, soil homogeneous
    Nco_0_0  = 2.483 * np.log(H_B) + 1.974  # angle = 0 deg 
    Nco_90_0 = 2.174 * np.log(H_B) + 3.391  # angle = 90 deg

    kBSh = k * B / Suh    # The degree of soil non-homogeneity

    f0  = np.where(H_B < 4, 1.77 * (H_B**0.3) - 1.289, 0.192 * H_B + 0.644)
    f90 = np.where(H_B < 4, 0.68 * (H_B**0.5) - 0.41 , 0.153 * H_B + 0.341)

    # non-homogeneity adjustment factor for anchor ultimate pullout capacity
    S_kB_0  = 1 - f0  * kBSh
    S_kB_90 = 1 - f90 * kBSh

    # Anchor Pullout capacity factor in weightless clay with breakaway base, soil nonhomogeneous
    Nco_0  = S_kB_0  * Nco_0_0  
    Nco_90 = S_kB_90 * Nco_90_0 

    # Anchor Pullout capacity factor in weightless clay with no breakaway base 
    Nco = Nco_0 + (Nco_90 - Nco_0) * (Bita/90)**2   


    # ----- ultimate anchor capacity factor -----

    # Ultimate Anchor capacity factor, soil homogeneous
    Nco_s_0_0  = np.where(2.90 * H_B + 6.02 <= 11.59, 2.90 * H_B + 6.02, 11.596)
    Nco_s_90_0 = np.where(2.72 * H_B + 4.02 <= 11.59, 2.72 * H_B + 4.02, 11.596)

    # non-homogeneity factor for anchor ultimate pullout capacity
    S_s_kB_0 = np.where(H_B <= 2, 1 + (0.8 - 0.3 * H_B) * kBSh - (0.383 * kBSh**1.36), 1)  # Angle = 0

    f90s = np.where(H_B <= 3, 0.267 * H_B, 0.6)
    S_s_kB_90 = 1 - f90s * kBSh  # Angle = 90

    # Anchor ultimate holding capacity in with breakaway base, soil nonhomogeneous
    Nco_s_0  = S_s_kB_0  * Nco_s_0_0
    Nco_s_90 = S_s_kB_90 * Nco_s_90_0

    # Anchor ultimate holding capacity in with no breakaway base, soil nonhomogeneous
    Nco_s = Nco_s_90 + (Nco_s_0 - Nco_s_90) * ((90-Bita) /90)**2


    # ----- final results -----

    Nc_final = np.minimum(Nco + (gamma * Hs)/Suh, Nco_s) # anchor pullout capacity factor [kN ??]

    qu = Nc_final * Suh  # The bearing pressure capacity of the anchor plate

    T = qu * (1-Los) * A  # The bearing tension force capacity of the anchor plate


    results = {}
    results['capacity'] = T  # capacity at specified loading angle
    results['vol'] = V
    results['L'] = B
    results['D'] = B
    results['t'] = t
    
    return results
