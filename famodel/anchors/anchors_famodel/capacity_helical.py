
import numpy as np

def getCapacityHelical(D, L, d, soil_type, gamma, alpha_star, Su0=None, k=None, phi=None):
    
    '''Calculate the inclined vertical load capacity of a helical pile in clay.
    The calculation is based on the soil properties and anchor geometry.  
    
    Parameters
    ----------
    D : float 
        Helix diameter. [m]
    L : float 
        Length shaft [m]
    d : float 
        Pile shaft diameter. [m]
    soil_type : string
        Select soil condition, 'clay' or 'sand' 
    gamma: float 
        Effective unit weight of the soil. [kN/m3] 
    Su0 : float 
        Undrained shear strength at the mudline [kPa]
    k : float 
        Undrained shear strength gradient [kPa/m]
    
    
    Returns
    -------
    Qu: float
       Maximum vertical capacity [kN]    
    '''
    
    # ----- Clay case -----
    if soil_type == 'clay':
        
        Nc = 6.0*(1 + 0.2*d/D); 
        Nc = np.where(Nc < 9, Nc, 9)
        # Su is calculated, at the depth of the helix minus one helical plate diameter
        # A reduction of 25% is applied for a moderately sensitive clay
        Qh = ((np.pi/4)*(D**2 - d**2)*Nc*(Su0 + k*(L - D)) + gamma*D)*0.75 
        Qs = np.pi*d*L*alpha_star*(Su0 + k*(L - D))
        Qu = Qh + Qs
                
    # ----- Sand case -----
    else:
    
        Nq = 0.5*(12*phi)**(phi/54) 
        Qh = (np.pi/4)*(D**2 - d**2)*Nq*gamma*L
        Qs = np.pi*d*L*alpha_star*gamma*L
        Qu = Qh + Qs
            
        
    resultsHelical = {}
    resultsHelical['Capacity'] = Qu   # Vertical capacity 

    return resultsHelical