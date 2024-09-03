
import numpy as np

def getCapacityPlate(A, beta, zlug, soil_type, gamma, Su0=None, k=None):
    
    '''Calculate the inclined load capacity of a plate in clay at a given depth.
    The calculation is based on the soil properties, anchor geometry and the angle of inclined load.
    The plate is assumed to be inclined perpendicular to the tension at the main padeye depth.     
    
    Parameters
    ----------
    A : float 
        Plate area, assumed to be square so that width B = sqrt(A). [m^2]
    beta : float
        Angle of the plate after keying process [deg]
    zlug: float
        Embedded depth of the main padeye [m]
    soil_type : string
        Specify 'sand' or 'clay'. This affects what other soil parameters are used.
    gamma: float 
        Effective unit weight of the soil [kN/m3] 
    Su0 : float 
        Undrained shear strength at the mudline [kPa]
    k : float 
        Undrained shear strength gradient [kPa/m]
   
    Returns
    -------
    Tmax: float
       Maximum capacity [kN] 
    '''
      

    Los=0.05                   # Key lost fraction due to the keying process, default 0.05 [-]
    rhos= 7850.0               # Steel density in air [kg/m3]
 
    B = round(np.sqrt(A),2)    # Anchor width (and length, approximated as square for now) [m]
    zlug_B = zlug/B               # Anchor depth range/ The width of the Plate 
    B_t = 40                   # Aspect ratio plate width to thickness, default is 40
    t = round(B/B_t, 2)        # Thickness of the plate, which it depends on the width [m]
    
    #t=np.sqrt(0.006*A)/4  
    V = round(A*t,2)           # Steel volume [m3]
    W = V*rhos                 # Plate weight [kg]
    Su = Su0 + k*zlug          # Undrained shear strength at plate depth

    # ----- anchor pullout capacity -----

    # Anchor Pullout capacity factor in weightless clay with breakaway base, soil homogeneous
    Nco_0_0  = 2.483*np.log(zlug_B) + 1.974  # angle = 0 deg 
    Nco_90_0 = 2.174*np.log(zlug_B) + 3.391  # angle = 90 deg

    kBSh = k*B/Su              # Degree of soil non-homogeneity

    f0  = np.where(zlug_B < 4, 1.77*(zlug_B**0.3) - 1.289, 0.192*zlug_B + 0.644)
    f90 = np.where(zlug_B < 4, 0.68*(zlug_B**0.5) - 0.41 , 0.153*zlug_B + 0.341)

    # Non-homogeneity adjustment factor for anchor ultimate pullout capacity
    S_kB_0  = 1 - f0 *kBSh
    S_kB_90 = 1 - f90*kBSh

    # Anchor Pullout capacity factor in weightless clay with breakaway base, soil nonhomogeneous
    Nco_0  = S_kB_0*Nco_0_0  
    Nco_90 = S_kB_90*Nco_90_0 

    # Anchor pullout capacity factor in weightless clay with no breakaway base 
    Nco = Nco_0 + (Nco_90 - Nco_0)*(beta/90)**2   

    # Uplift bearing capacity factor, soil homogeneous
    Nco_s_0_0  = np.where(2.90*zlug_B + 6.02 <= 11.59, 2.90*zlug_B + 6.02, 11.596)
    Nco_s_90_0 = np.where(2.72*zlug_B + 4.02 <= 11.59, 2.72*zlug_B + 4.02, 11.596)

    # ----- ultimate anchor capacity factor -----
    
    # Non-homogeneity factor for anchor ultimate pullout capacity
    S_s_kB_0 = np.where(zlug_B <= 2, 1 + (0.8 - 0.3*zlug_B)*kBSh - (0.383*kBSh**1.36), 1)  # Angle = 0

    f90s = np.where(zlug_B <= 3, 0.267*zlug_B, 0.6)
    S_s_kB_90 = 1 - f90s*kBSh  # Angle = 90

    # Anchor ultimate holding capacity in with breakaway base, soil nonhomogeneous
    Nco_s_0  = S_s_kB_0 *Nco_s_0_0
    Nco_s_90 = S_s_kB_90*Nco_s_90_0

    # Anchor ultimate holding capacity in with no breakaway base, soil nonhomogeneous
    Nco_s = Nco_s_90 + (Nco_s_0 - Nco_s_90)*((90 - beta)/90)**2


    # ----- final results -----
    Nc_final = np.minimum(Nco + (gamma*zlug)/Su, Nco_s) # anchor pullout capacity factor [kN]
    qu = Nc_final*Su  # The bearing pressure capacity of the anchor plate
    Tmax = round(qu*(1 - Los)*A,2)  # The bearing tension force capacity of the anchor plate
    Hmax = Tmax*np.cos(np.deg2rad(beta))
    Vmax = Tmax*np.sin(np.deg2rad(beta))

    resultsPlate = {}
    resultsPlate['Capacity'] = Tmax                 # Capacity at specified loading angle
    resultsPlate['Horizontal max.'] = Hmax          # Maximum horizontal capacity in sand
    resultsPlate['Vertical max.'] = Vmax            # Maximum vertical capacity in sand
    
    return resultsPlate