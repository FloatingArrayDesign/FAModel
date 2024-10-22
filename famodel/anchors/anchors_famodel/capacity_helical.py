
import numpy as np

def getCapacityHelical(D, L, d, soil_type, gamma, Su0=None, k=None, phi=None, Dr=None):
    
    '''Calculate the inclined vertical load capacity of a helical pile in clay.
    The calculation is based on the soil properties and anchor geometry.  
    
    Parameters
    ----------
    D : float 
        Helix diameter [m]
    L : float 
        Length shaft [m]
    d : float 
        Pile shaft diameter [m]
    soil_type : string
        Select soil condition, 'clay' or 'sand' 
    gamma: float 
        Effective unit weight of the soil [kN/m3] 
    Su0 : float 
        Undrained shear strength at the mudline (clay only) [kPa]
    k : float 
        Undrained shear strength gradient (clay only) [kPa/m]
    phi : float
        Angle of internal friction (sand only) [deg]
    Dr : float
        Relative density of the soil (%) (sand only) [-] 
    
    
    Returns
    -------
    Qu: float
       Maximum vertical capacity [kN]    
    '''
    
    rhos= 78.50                      # Dry steel unit weight (kN/m3)
    t = (6.35 + D*20)/1e3            # Suction pile wall thickness (m), API RP2A-WSD
    
    # Dry and wet mass of the pile    
    def PileWeight(Len, Dia1, Dia2, tw, rho):
        Wp = ((np.pi/4)*((Dia1**2 - (Dia1 - 2*tw)**2)*Len + (np.pi/4)*Dia2**2*tw))*rho
        return Wp
    # Define alpha coefficient (clay)
    sigma_v_eff = gamma*zlug         # Effective soil stress (kN/m2)
    psi_val = Su_av_L/sigma_v_eff    # Su/p0' for point in question (API DP 2A-WSD)   
    if psi_val <= 1.0:
        alpha = min(0.5*psi_val**-0.50, 1)
    else:
        alpha = min(0.5*psi_val**-0.25, 1)
        
    # Define delta as a function of Dr (sand)
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
        
    Wp = PileWeight(L, D, d, t, rhos) 
    
    # ----- Clay case -----
    if soil_type == 'clay':       
        Nc = 6.0*(1 + 0.2*d/D); 
        Nc = np.where(Nc < 9, Nc, 9)
        # Su is calculated, at the depth of the helix minus one helical plate diameter
        # A reduction of 25% is applied for a moderately sensitive clay
        Qh = ((np.pi/4)*(D**2 - d**2)*Nc*(Su0 + k*(L - D)) + gamma*D)*0.75 
        Qs = np.pi*d*L*alpha*(Su0 + k*(L - D))
        Qu = Qh + Qs
                
    # ----- Sand case -----
    else:   
        delta = calc_delta(Dr)
        Nq = 0.5*(12*phi)**(phi/54) 
        Qh = (np.pi/4)*(D**2 - d**2)*Nq*gamma*L
        Qs = np.pi*d*L*delta*gamma*L
        Qu = Qh + Qs
                  
    resultsHelical = {}
    resultsHelical['Capacity'] = Qu          # Vertical capacity 
    resultsHelical['Pile weight'] = Wp       # Dry weight of the helical pile (kN)

    return resultsHelical