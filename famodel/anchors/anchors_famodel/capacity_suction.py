
import yaml      
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
#from famodel.anchors.capacity_load import getAnchorLoad
from famodel.anchors.capacity_load import getAnchorLoadDNV

def getCapacitySuction(D, L, zlug, H, V, soil_type, gamma, Su0=None, k=None, phi=None, Dr=None, plot=True):
    
    '''Calculate the inclined load capacity of a suction pile in sand or clay following S. Kay methodology.
    The calculation is based on the soil properties, anchor geometry and inclined load.  
    
    Parameters
    ----------
    D : float 
        Suction pile diameter [m]
    L : float 
        Suction anchor length [m]
    zlug: float
        Embedded depth of the main padeye [m]
    soil_type : string
        Select soil condition, 'clay' or 'sand'              
    gamma: float 
        The effective unit weight of the soil. [kN/m3]
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
    Hmax : float 
        Maximum horizontal capacity [kN]
    Vmax : float 
        Maximum vertical capacity [kN]
    '''  
            
    lambdap = L/D; m = -2/3;         # Suction pile slenderness ratio
    t = (6.35 + D*20)/1e3            # Suction pile wall thickness (m), API RP2A-WSD
    rlug = D/2                       # Radial position of the lug
    thetalug = 5                     # Angle of tilt misaligment, default is 5. (deg)
    psilug = 7.5                     # Angle of twist misaligment, default is 7.5. (deg)
    rhows = 66.90                    # Submerged steel specific weight (kN/m3)
    rhow = 10                        # Water specific weight (kN/m3) 
    
    # Outer and inner surface of the pile skirt
    def PileSurface(Len, Dia):
        Sp = np.pi*Dia*Len
        return Sp    
    # Dry and wet mass of the pile    
    def PileWeight(Len, Dia, tw, rho):
        Wp = ((np.pi/4)*((Dia**2 - (Dia - 2*tw)**2)*Len + (np.pi/4)*Dia**2*tw))*rho
        return Wp 
    # Mass of the soil plug      
    def SoilWeight(Len, Dia, tw, gamma_soil): 
        Wsoil =(np.pi/4)*(Dia - 2*tw)**2*Len*gamma_soil
        return Wsoil
    # Tilt and twist effects due to installation misaligments
    def rlugTilt(r, z, theta):
        R = r*np.cos(np.deg2rad(theta)) - z*np.sin(np.deg2rad(theta))
        return R   
    def zlugTilt(r, z, theta):
        Z = r*np.sin(np.deg2rad(theta)) + z*np.cos(np.deg2rad(theta))
        return Z
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
       
    if soil_type == 'clay':
        # Definitions for cohesive soils
        Nc = min (6.2*(1 + 0.34*np.arctan(lambdap)),9)   # End-bearing capacity factor
        ez = (Su0*L**2/2 + k*L**3/3)/(Su0*L + k*L**2/2)
        Np_fixed = 10.25; Np_free = 4                    # From Np vs L/D chart from CAISSON_VHM
        Su_av_L = Su0 + k*zlug                           # Undrained shear strength values (average) 
        Su_tip = Su0 + k*L                               # Undrained shear strength values (tip)
        sigma_v_eff = gamma*zlug                         # Effective soil stress (kN/m2)
        psi_val = Su_av_L/sigma_v_eff                    # Su/p0' for point in question (API DP 2A-WSD)
        #zlug = ez                                       # Optimized depth of the lug 
 
        if psi_val <= 1.0:
            alpha = min(0.5*psi_val**-0.50, 1)
        else:
            alpha = min(0.5*psi_val**-0.25, 1)
 
        Hmax = Np_fixed*L*D*Su_av_L; 
        H0 = Np_free*L*D*Su_av_L;
        Mmax = Np_fixed*L*L*D*Su_av_L; 
        
        # M modifies the Hmax capacity
        M = - V*rlugTilt(rlug,zlug,thetalug) - H*(zlugTilt(rlug,zlug,thetalug) - ez)
        def f(Hmax):
             return m*(Hmax/(L*D*(Su0 + k*zlug)) - Np_fixed) + M*(Hmax/(L*D*(Su0 + k*zlug))/(Hmax*L))
        Hmax = fsolve(f,5);
        
        # Torsion capacity
        Fo = PileSurface(L, D)*alpha*Su_av_L
        To = Fo
        Ti = PileSurface(L,(D - 2*t))*alpha*Su_av_L
        Tbase = np.pi*D**3*Su_tip/12
        Tmax = min(To + Ti, To + Tbase) 
        
        # Introduce twist effects due to installation misaligment
        T = H*rlug*np.sin(np.deg2rad(psilug))
        nhuT = T/Tmax; nhuV = H/Fo;
        nhuVstar = np.sqrt(nhuV**2 - nhuT**2)
        alphastar = alpha*(nhuVstar/nhuV)
        
        # "Plugged" (Reverse end bearing capacity - passive suction) 
        Vmax1 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alphastar*Su_av_L + Nc*Su_tip*(np.pi/4)*D**2)
        # "Coring"        
        Vmax2 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alphastar*Su_av_L + PileSurface(L,(D - 2*t))*alphastar*Su_av_L)
        # "Leaking"        
        Vmax3 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alphastar*Su_av_L + SoilWeight(L, D, t, gamma))                  
        Vmax = min(Vmax1, Vmax2, Vmax3)
               
    elif soil_type == 'sand':
        # Definition for non-cohesive soils
        Nq = np.e**(np.pi*np.tan(np.radians(phi)))*np.tan(np.radians(45) + np.radians(phi)/2)**2 # Lateral-bearing capacity factor
        sigma_av_L = gamma*L/2                        # Effective stress (average)
        sigma_tip = gamma*L                           # Effective stress (tip)
        Hmax = 0.5*D*Nq*gamma*L**2

        M = - V*rlugTilt(rlug,zlug,thetalug) - H*(zlugTilt(rlug,zlug,thetalug) - zlug)
        
        # Torsion capacity
        delta = calc_delta(Dr)
        To = PileSurface(L, D)*delta*sigma_av_L
        Ti = PileSurface(L, (D -2*t))*delta*sigma_av_L
        Tbase = np.pi*D**3*sigma_tip/12
        Tmax = min(To + Ti, To + Tbase) 
        
        # Introduce twist effects due to installation misaligment
        T = H*rlug*np.sin(np.deg2rad(psilug))
        Fo = delta*sigma_av_L*L*np.pi*D
        nhuT = T/Tmax; nhuV = H/Fo;
        nhuVstar = np.sqrt(nhuV**2 - nhuT**2)
        deltastar = delta*(nhuVstar/nhuV)
    
        # "Coring"        
        Vmax2 = PileWeight(L, D, t, rhows) + PileSurface(L, D)*deltastar*sigma_av_L + PileSurface(L,(D - 2*t))*deltastar*sigma_av_L
        # "Leaking"        
        Vmax3 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*deltastar*sigma_av_L + SoilWeight(L, D, t, gamma))
        Vmax = min(Vmax2, Vmax3)
        # def y(depth):
            # return np.e**(-depth) - 1 + depth
        # Ze = D/(4*7); Zi = D/(4*5)    
        # Vmax = 7*gamma*Ze**2*y(L/Ze)*PileSurface(L, D)/L + 5*gamma*Zi**2*y(L/Zi)*PileSurface(L,(D - 2*t))/L
     
    # Pile weight (inc. stiffening plus vent) assessed as a factor
    Wp = 1.10*PileWeight(L, D, t, (rhows + rhow)) 
    # Submerged weight of the soil plug
    Ws = SoilWeight(L, D, t, gamma)
    
    # Capacity envelope
    aVH = 0.5 + lambdap; bVH = 4.5 + lambdap/3 
    # print('Env. exp = ' +str(aVH)+'   '+str(bVH))
    UC = (H/Hmax)**aVH + (V/Vmax)**bVH      
    x = np.cos(np.linspace (0, np.pi/2, 100))
    y = (1 - x**bVH)**(1/aVH)
    X = Hmax*x; Y = Vmax*y
    if plot:
        plt.plot(X, Y, color = 'b')
        plt.plot(H, V, color = 'r')
            
        # Set labels and title
        plt.xlabel('Horizontal capacity [kN]')
        plt.ylabel('Vertical capacity [kN]')
        plt.suptitle('VH suction pile capacity envelope')
        plt.axis([0, 1.3*max(X[0], H), 0, 1.3*max(Y[-1], V)]) 
        plt.grid(True)
        plt.show()
    
    resultsSuction = {}
    if soil_type == 'clay':
        resultsSuction['Horizontal max.'] = Hmax[0] # Maximum horizontal capacity in clay
    elif soil_type == 'sand':   
        resultsSuction['Horizontal max.'] = Hmax    # Maximum horizontal capacity in sand 
    resultsSuction['Vertical max.'] = Vmax          # Maximum vertical capacity 
    if soil_type == 'clay':
        resultsSuction['UC'] = UC[0]                # Unity check in clay
    elif soil_type == 'sand':
        resultsSuction['UC'] = UC                   # Unity check in sand
    resultsSuction['Weight Pile'] = Wp              # Dry weight of the suction pile (kN)
    resultsSuction['Weight Soil'] = Ws              # Submerged weight of the soil plug (kN)
    resultsSuction['t'] = t                         # Pile thikness in [m]
    
    return resultsSuction
        
def getCapacitySuctionSimp(D, L, zlug, H, V, gamma, Su0, k, alpha):

    '''
    Parameters
    ----------
    D : float 
        Suction pile diameter [m]
    L : float 
        Suction anchor length [m]
    Tm : float 
        Mooring line load at mudlevel [kN]
    thetam : float 
        Mooring line angle at mudlevel [deg]
    zlug : float 
        Embedded depth of the lug [m]
    gamma: float
    
    Su0 : float 
        Undrained shear strength at the mudline (clay only)[kPa]
    k : float 
        Undrained shear strength gradient (clay only) [kPa/m]
    alpha : float 
        Skin friction coefficient (outer and inner - clay only) [-] 
    rhows : float
        Submerged steel density [t/m3]

    Returns
    -------
    Hmax : float 
        Maximum horizontal capacity [kN]
    Vmax : float 
        Maximum vertical capacity [kN]
    UC: float
        Capacity unity check for given load [-]
    ''' 
    
    lambdap = L/D;                       # Suction pile slenderness ratio   
    t = 10*D/1e3                         # Thickness of the pile
    Np_fixed = 10.25;                    # From Np vs L/D chart from CAISSON_VHM
    rhows=66.90
    
    Su_av_L = Su0 + k*zlug;              # Undrained shear strength values (average)
    Su_tip = Su0 + k*L;                  # Undrained shear strength values (tip)
    Nc = min (6.2*(1 + 0.34*np.arctan(lambdap)),9)  # End-bearing capacity factor
    
    # Outer and inner surface of the pile skirt
    def PileSurface(Len, Dia):
        Sp = np.pi*Dia*Len
        return Sp    
    # Dry and wet mass of the pile    
    def PileWeight(Len, Dia, tw, rho):
        Wp = ((np.pi/4)*((Dia**2 - (Dia - 2*tw)**2)*Len + (np.pi/4)*Dia**2*tw))*rho
        return Wp 
    # Weight of the soil plug      
    def SoilWeight(Len, Dia, tw, gamma_soil): 
        Wsoil =(np.pi/4)*(Dia - 2*tw)**2*Len*gamma_soil
        return Wsoil
    
    Hmax = Np_fixed*L*D*Su_av_L
    # "Plugged" (Reverse end bearing capacity - passive suction) 
    Vmax1 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alpha*Su_av_L + Nc*Su_tip*(np.pi/4)*D**2)
    # "Coring"        
    Vmax2 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alpha*Su_av_L + PileSurface(L,(D - 2*t))*alpha*Su_av_L)
    # "Leaking"        
    Vmax3 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alpha*Su_av_L + SoilWeight(L, D, t, gamma)) 
    # print(Vmax1); print(Vmax2); print(Vmax3)                  
    Vmax = min(Vmax1,Vmax2,Vmax3)
    
    # Submerged pile weight (inc. stiffening plus vent) assessed as a factor
    Wp = 1.00*PileWeight(L, D, t, (rhows)) 
    # Submerged weight of the soil plug
    Ws = SoilWeight(L, D, t, gamma) 
    
    # H = Tm*np.cos(np.deg2rad(thetam)); V = Tm*np.sin(np.deg2rad(thetam))
    
    # Capacity envelope
    aVH = 0.5 + lambdap; bVH = 4.5 + lambdap/3 
    # print('Env. exp =' +str(aVH)+str(bVH))
    UC = (H/Hmax)**aVH + (V/Vmax)**bVH

    x = np.cos(np.linspace (0,np.pi/2,1000))
    y = (1 - x**bVH)**(1/aVH)
    X = Hmax*x; Y = Vmax*y

    #H_good = Hmax*np.exp(np.log(0.5)/aVH)
    #V_good = Vmax*np.exp(np.log(0.5)/bVH)
   
    resultsSuctionSimp = {}
    resultsSuctionSimp['Horizontal max.'] = Hmax    # Capacity at specified loading angle
    resultsSuctionSimp['Vertical max.'] = Vmax      # Capacity at specified loading angle
    resultsSuctionSimp['UC'] = UC                   # Unity check
    resultsSuctionSimp['Weight Pile'] = Wp          # in kN
    resultsSuctionSimp['Weight Soil'] = Ws          # in kN
    resultsSuctionSimp['t'] = t
    
    return resultsSuctionSimp    