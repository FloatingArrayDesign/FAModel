"""Screw pile capacity calculation functions in sandy/clayey soils, designed
especially for computing the uplift load perpendicular to the anchor plate.
Lead author: Felipe Moreno. 
"""

import numpy as np

def getCapacityScrew(D, L, d, Su0=1.4, k=1.61, soil_type='clay', 
                     gamma=7.3, phi=38, alpha_star=0.60, alphac=0.35):
    '''Calculate the inclined load capacity of a suction-embedded plate anchor in clay.
    The calculation is based on the soil properties, anchor geometry, and the angle of inclined load.  
    
    Parameters
    ----------
    D : float 
        Helix diameter. [m]
    L : float 
        Length shaft [m]
    d : float 
        Pile shaft diameter. [m]
    S : float
        Distance between helices. [deg]
    p : float
        Helix pith distance. [m]
    n : float
        Number of helices, default 1. [-]
    soil_type : string
        Specify 'sand' or 'clay'. This affects what other soil parameters are used.
    gamma: float 
        Effective unit weight of the sand soil, default is 9 for sand, 4.7 for clay. [kN/m3] 
    phip : float 
        Sand peak friction (sand only), default is 39. [kPa]
    psip : float
        Dilatancy angle (sand only), default 9. [kPa/m]
    
    
    Returns
    -------
    
    '''
    
    # ----- Clay case -----
    if soil_type == 'clay':
        
        Nc = 6.0*(1 + 0.2*d/D); Nc = np.where(Nc < 9, Nc, 9); print(Nc)
        # Su is calculated, at the depth of the helix minus one helical plate diameter
        # A reduction of 25% is applied for a moderately sensitive clay
        Qh = ((np.pi/4)*(D**2 - d**2)*Nc*(Su0 + k*(L - D)) + gamma*D) 
        Qs = np.pi*d*L*alpha_star*(Su0 + k*(L - D))
        Qu = Qh + Qs
                
    # ----- Sand case -----
    else:
    
        Nq = 0.5*(12*phi)**(phi/54) 
        Qh = (np.pi/4)*(D**2 - d**2)*Nq*gamma*L
        Qs = np.pi*d*L*alpha_star*gamma*L
        Qu = Qh + Qs
        
    # #Calculation of the uplift factors
    # Fps = np.tan(np.deg2rad(phip)) + np.cos(np.deg2rad(phip) - np.deg2rad(psip))*np.tan(np.deg2rad(phip - psip)) 
    # Fs1 = 2*Fps
    # Fs2 = 4/3*Fps*np.tan(np.deg2rad(psip))
    
    
    # #----- Anchor pullout capacity -----
    # Fu = (1 + Fs1*(L/D) + Fs2*(L/D)**2)*gamma*np.pi*D**2*L/4                 
        
    results = {}
    results['capacity1'] = Qu  # capacity at specified loading angle
    #results['capacity2'] = Fu  # capacity at specified loading angle

    
    return results


if __name__ == '__main__':
 
    
    ''' Testing the function in one case of the screw = 10 m2, 
    the aspect ratio of the plate width and its thickness, default is 40, all other parameters are default.'''
 
    results = getCapacityScrew(5.01, 25.1, 1.00)
    
    print('********************* One Case Test Result********************')

    print('Load capacity,     '          , results['capacity1'], '[kN]')
    #print('Load capacity,     '          , results['capacity2'], '[kN]')

    print('**************************************************************') 