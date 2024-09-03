import numpy as np

def getCapacityPlate(A, Hs=20, soil_type='clay', gamma=0, Su0=2.39, k=1.41, beta=30):
    
    '''Calculate the inclined load capacity of a suction-embedded plate anchor in clay.
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
    beta : float
        The Angle of the plate after Keying Process, default 30. [deg]
    Hs   : float
        The Specified of the plate embedded depth, default 20. [m]
    Los  : float
        The undrained shear strength lost fraction due to the keying process, default 0.05.
    
    Returns
    -------
    
    '''
    
    if soil_type == 'sand':
        raise Exception('Only clay soil is supported so far.')
    

    Los=0.05                   # Key loss factor
    rhos= 7850.0               # The steel density in kg/m3
 
    B = round(np.sqrt(A),2)    # Anchor width (and length, approximated as square for now) [m]
    H_B = Hs/B                 # Anchor depth range/ The width of the Plate 
    B_t = 40                   # Aspect ratio plate width to thickness
    t = round(B/B_t,2)         # Thickness of the plate, which it depends on the width [m]
    
    #t=np.sqrt(0.006*A)/4  
    V = round(A*t,2)           # Steel volume [m3]
    W = V*rhos                 # Plate weight [kg]
    Su = Su0 + k*Hs            # Undrained shear strength at plate depth


    # ----- anchor pullout capacity -----

    # Anchor Pullout capacity factor in weightless clay with breakaway base, soil homogeneous
    Nco_0_0  = 2.483*np.log(H_B) + 1.974  # angle = 0 deg 
    Nco_90_0 = 2.174*np.log(H_B) + 3.391  # angle = 90 deg

    kBSh = k*B/Su              # Degree of soil non-homogeneity

    f0  = np.where(H_B < 4, 1.77*(H_B**0.3) - 1.289, 0.192*H_B + 0.644)
    f90 = np.where(H_B < 4, 0.68*(H_B**0.5) - 0.41 , 0.153*H_B + 0.341)

    # Non-homogeneity adjustment factor for anchor ultimate pullout capacity
    S_kB_0  = 1 - f0 *kBSh
    S_kB_90 = 1 - f90*kBSh

    # Anchor Pullout capacity factor in weightless clay with breakaway base, soil nonhomogeneous
    Nco_0  = S_kB_0*Nco_0_0  
    Nco_90 = S_kB_90*Nco_90_0 

    # Anchor pullout capacity factor in weightless clay with no breakaway base 
    Nco = Nco_0 + (Nco_90 - Nco_0)*(Beta/90)**2   

    # Uplift bearing capacity factor, soil homogeneous
    Nco_s_0_0  = np.where(2.90*H_B + 6.02 <= 11.59, 2.90*H_B + 6.02, 11.596)
    Nco_s_90_0 = np.where(2.72*H_B + 4.02 <= 11.59, 2.72*H_B + 4.02, 11.596)

    # ----- ultimate anchor capacity factor -----
    
    # Non-homogeneity factor for anchor ultimate pullout capacity
    S_s_kB_0 = np.where(H_B <= 2, 1 + (0.8 - 0.3*H_B)*kBSh - (0.383*kBSh**1.36), 1)  # Angle = 0

    f90s = np.where(H_B <= 3, 0.267*H_B, 0.6)
    S_s_kB_90 = 1 - f90s * kBSh  # Angle = 90

    # Anchor ultimate holding capacity in with breakaway base, soil nonhomogeneous
    Nco_s_0  = S_s_kB_0 *Nco_s_0_0
    Nco_s_90 = S_s_kB_90*Nco_s_90_0

    # Anchor ultimate holding capacity in with no breakaway base, soil nonhomogeneous
    Nco_s = Nco_s_90 + (Nco_s_0 - Nco_s_90)*((90 - beta)/90)**2


    # ----- final results -----

    Nc_final = np.minimum(Nco + (gamma*Hs)/Su, Nco_s) # anchor pullout capacity factor [kN]
    qu = Nc_final*Su  # The bearing pressure capacity of the anchor plate
    T = round(qu*(1 - Los)*A,2)  # The bearing tension force capacity of the anchor plate


    resultsPlate = {}
    resultsPlate['capacity'] = T  # capacity at specified loading angle
    resultsPlate['vol'] = V
    resultsPlate['L'] = B
    resultsPlate['D'] = B
    resultsPlate['t'] = t
    
    return results


if __name__ == '__main__':
 
    
    ''' Testing the function in one case of the plate with area = 10 m2, 
    the aspect ratio of the plate width and its thickness, default is 40, all other parameters are default.'''
 
    results = getCapacityPlate(10)
    print('********************* One Case Test Result********************')

    print('Anchor length,              ' , resultsPlate['L'], '[m]') 
    print('Anchor diameter,            ' , resultsPlate['D'], '[m]') 
    print('Anchor thickness,           ' , resultsPlate['t'], '[m]')
    print('Anchor steel volume,        ' , resultsPlate['vol'], '[m3]') 
    print('Inclined load capacity,     ' , resultsPlate['capacity'], '[kN]') 

    print('**************************************************************') 
