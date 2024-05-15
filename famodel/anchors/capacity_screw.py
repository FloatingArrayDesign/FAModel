"""Screw pile capacity calculation functions in sandy/clayey soils, designed
especially for computing the uplift load perpendicular to the anchor plate.
Lead author: Felipe Moreno. 
"""

import numpy as np

def getCapacityScrew(D,L,d,Su0,dSu,soil_type='sand', 
                     gamma=10,phi=50,alphas=0.35,alphac=0.35):
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
        
        Nc = 6.0*(1 + 0.2*d/D); Nc = np.where(Nc < 9, Nc, 9)
        Qh = (np.pi/4)*(D**2 - d**2)*(Su0 + dSu/3) +  gamma*D  
        Qs = np.pi*d*L*alphas*(Su0 + dSu/3)
        Qu = Qh + Qs
                
    # ----- Sand case -----
    else:
    
        Nq = 0.5*(12*phi)**(phi/54) 
        Qh = (np.pi/4)*(D**2 - d**2)*gamma*L*Nq
        Qs = np.pi*d*L*alphas*gamma*L
        Qu = Qh + Qs
        
        # Calculation of the uplift factors
        # Fps = np.tan(np.deg2rad(phip)) + np.cos(np.deg2rad(phip) - np.deg2rad(psip))*np.tan(np.deg2rad(phip - psip)) 
        # Fs1 = 2*Fps
        # Fs2 = 4/3*Fps*np.tan(np.deg2rad(psip))
        
        # ----- Anchor pullout capacity -----
        # Fu = (1 + Fs1*(L/D) + Fs2*(L/D)**2)*gamma*np.pi*D**2*L/4                 
        
    results = {}
    results['capacity'] = Qu  # capacity at specified loading angle

    
    return results


if __name__ == '__main__':
 
    
    ''' Testing the function in one case of the screw = 10 m2, 
    the aspect ratio of the plate width and its thickness, default is 40, all other parameters are default.'''
 
    results = getCapacityScrew(0.5, 4.3)
    print('********************* One Case Test Result********************')

    print('Load capacity,     '          , results['capacity'], '[kN]') 

    print('**************************************************************') 