
import yaml     
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from matplotlib import cm
from mpl_toolkits import mplot3d 
   
def getCapacityTorpedo(D1, D2, L1, L2, zlug, soil_type, Su0, k, alpha):
    
    '''Calculate the inclined load capacity of a torpedo pile in clay following S. Kay methodology.
    The calculation is based on the holding capacity of the suction pile as if it fully was embedded in soil.  
    
    Parameters
    ----------
    D1 : float 
        Torpedo pile wing diameter. [m]
    D2 : float    
        Torpedo pile shaft diameter. [m]
    L1 : float 
        Torpedo pile wing length. [m]
    L2 : float 
        Torpedo pile shaft length, excluding wing. [m]
    zlug : float 
        Torpedo pile embedded depth at main padeye elevation. [m]        
    soil_type : string
        Select soil condition, 'clay' or 'sand' 
    gamma: float 
        The effective unit weight of the soil. [kN/m3] 
    Su0 : float 
        The Undrained shear strength at the mudline. [kPa]
    k : float 
        The Undrained shear strength gradient. [kPa/m]
    alpha : float
        The skin friction coefficient. [-]   
    
    Returns
    -------
    Hmax : float 
        Maximum horizontal capacity [kN]
    Vmax : float 
        Maximum vertical capacity [kN]
    '''

    #m = -5;
    L = L1 + L2;
    Dstar = (D1*L1 + (D1 + 2*D2)*L2)/L              # Plane 1 (four fins)
    #Dstar = (D1*L1 + np.sqrt(2)*(D1/2 + D2)*L2)/L  # Plane 2 (four fins)    
    #rlug = D2/2; zlug = zlug; 
    lambdap = L/Dstar; print('lambdap = ' +str(lambdap))
    a = zlug; b = zlug + L1; c = zlug + L1 + L2;
    Wp = 850                                        # Weight of the pile. [kN]
    
    # Dry and wet mass of the pile
    def PileSurface(Len1, Len2, Dia1, Dia2):
        Sp = np.pi*Dia1*(Len1 + Len2) + 8*Len2*Dia2*0.9
        return(Sp)  
    # Dry and wet mass of the pile    
    def PileWeight(Len1, Len2, Dia1, Dia2, tw, rho):
        Wp = ((np.pi/4)*(Dia1**2 - (Dia1 - 2*tw)**2)*(Len1 + Len2) + 4*Len2*Dia2*tw)*rho
        return(Wp)
    
    ez_Su_den = D1*Su0*(b - a) + 0.5*D1*k*(b**2 - a**2) + D2*Su0*(c - b) + 0.5*D2*k*(c**2 - b**2)
    ez_Su_num = D1*Su0*(a**2 - a*b) + 0.33*D1*k*(b**3 - a**3) + b**2*(0.5*D1*Su0 - 0.5*D1*a*k) - a**2*(0.5*D1*Su0 - 0.5*D1*a*k)\
        + D2*Su0*(b**2 - b*c) + 0.33*D2*k*(c**3 - b**3) + c**2*(0.5*D2*Su0 - 0.5*D2*b*k) - b**2*(0.5*D2*Su0 - 0.5*D2*b*k)
    ez_Su = ez_Su_num/ez_Su_den
    ez_Su_L = ez_Su/L
    print('ez_Su = ' +str(ez_Su))
    Np_free = 3.45     # From Np vs L/D chart from CAISSON_VHM
                 
    Hmax = L*Dstar*Np_free*(Su0 + k*(zlug + ez_Su))
    print('Hmax = ' +str(Hmax))
    Vmax = PileSurface(L1, L2, D1, D2)*alpha*(Su0 + k*(zlug + ez_Su)) + Wp
    print('Vmax = ' +str(Vmax))
       
    #aVH = 0.5 + L/Dstar; bVH = 4.5 - L/(3*Dstar)
    aVH = 4.5 + L/(2*Dstar); bVH = 3.5 - L/(4*Dstar) 
    #H = Ta*np.cos(np.deg2rad(resultsLoad['angle'])); V = Ta*np.sin(np.deg2rad(resultsLoad['angle']))
    #UC = (H/Hmax)**aVH + (V/Vmax)**bVH
    
    deg = [0, 15, 30, 45, 60, 75, 90]
    x = np.cos(np.deg2rad(deg))
    y = (1 - x**bVH)**(1/aVH)
    X = Hmax*x/1e3; Y = Vmax*y/1e3  # in MN
       
    resultsTorpedo = {}
    resultsTorpedo['Horizontal max.'] = Hmax #Hmax[0]    # Capacity at specified loading angle
    resultsTorpedo['Vertical max.'] = Vmax               # Capacity at specified loading angle
   
    return resultsTorpedo