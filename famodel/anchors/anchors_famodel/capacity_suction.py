# # -*- coding: utf-8 -*-
# """
# Created on Wed May 29 15:53:52 2024

# @author: fmoreno
# """

import yaml      # Allow access to config file for user inputs
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
#from famodel.anchors.capacity_load import getAnchorLoad
from famodel.anchors.capacity_load import getAnchorLoadDNV


def getCapacitySuction(D, L, soil_type='clay', gamma = 8.75, Su0=15.0, k=2.0, alpha=0.515, phi=30, beta=0.29, rhows=66.90):
    
    '''Calculate the inclined load capacity of a suction pile in sand or clay.
    The calculation is based on the soil properties, anchor geometry and inclined load.  
    
    Parameters
    ----------
    D : float 
        Suction pile diameter [m]
    L : float 
        Suction anchor length [m]
    thetalug : float 
        Angle of tilt misaligment, default is 5. [deg]
    psilug : float 
        Angle of twist misaligment, default is 7.5. [deg]
    safety_factor: float
        Specify 'yes' or 'no'. This affects load and resistance factors.
    soil_type : string
        Specify 'clay' or 'sand'. This affects what other soil parameters.               
    Su0 : float 
        Undrained shear strength at the mudline (clay only) [kPa]
    k : float 
        Undrained shear strength gradient (clay only) [kPa/m]
    alpha : float 
        Skin friction coefficient (outer and inner - clay only), default is 0.7 [-] 
    phi : float
        Angle of internal friction (sand only) [deg]
    beta : float
        Skin friction coefficient (sand only), default is 0.46 [-]
    rhows : float
        Submerged steel specific weight [kN/m3]
    
    Returns
    -------
    Hmax : float 
        Maximum horizontal capacity [kN]
    Vmax : float 
        Maximum vertical capacity [kN]
    UC: float
        Capacity unity check for given load [-]
    '''  
            
    lambdap = L/D; m = -2/3;         # Suction pile slenderness ratio
    t = 10*D/1e3                     # Thickness of the pile
    rlug = D/2          # Radial position of the lug
    
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
       
    if soil_type == 'clay':
        # Setting default gamma values per soil type [kN/m3]
        gamma = 8.75
        # Definitions for cohesive soils
        Nc = min (6.2*(1 + 0.34*np.arctan(lambdap)),9)   # End-bearing capacity factor
        ez = (Su0*L**2/2 + k*L**3/3)/(Su0*L + k*L**2/2)
        Np_fixed = 10.25; Np_free = 4                    # From Np vs L/D chart from CAISSON_VHM
        Su_av_L = Su0 + k*zlug; Su_tip = Su0 + k*L       # Undrained shear strength values (average, tip)
        #zlug = ez                                       # Optimized depth of the lug 
        
        Hmax = Np_fixed*L*D*Su_av_L; H0 = Np_free*L*D*Su_av_L;
        Mmax = Np_fixed*L*L*D*Su_av_L; 
        
        # M modifies the Hmax capacity
        M = - resultsLoad['V']*rlugTilt(rlug,zlug,thetalug) - resultsLoad['H']*(zlugTilt(rlug,zlug,thetalug) - ez)
        def f(Hmax):
             return m*(Hmax/(L*D*(Su0 + k*zlug)) - Np_fixed) + M*(Hmax/(L*D*(Su0 + k*zlug))/(Hmax*L))
        Hmax = fsolve(f,5);
        
        # Torsion capacity
        Fo = PileSurface(L, D)*alpha*Su_av_L
        To = Fo; print(To)
        Ti = PileSurface(L,(D - 2*t))*alpha*Su_av_L; print(Ti)
        Tbase = np.pi*D**3*Su_tip/12
        Tmax = min(Ti + To, To + Tbase) 
        
        # Introduce twist effects due to installation misaligment
        T = resultsLoad['H']*rlug*np.sin(np.deg2rad(psilug))
        nhuT = T/Tmax; nhuV = resultsLoad['H']/Fo;
        nhuVstar = np.sqrt(nhuV**2 - nhuT**2)
        alphastar = alpha*(nhuVstar/nhuV)
        
        # "Plugged" (Reverse end bearing capacity - passive suction) 
        Vmax1 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alphastar*Su_av_L + Nc*Su_tip*(np.pi/4)*D**2)/FoSsoil
        # "Coring"        
        Vmax2 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alphastar*Su_av_L + PileSurface(L,(D - 2*t))*alphastar*Su_av_L)/FoSsoil
        # "Leaking"        
        Vmax3 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alphastar*Su_av_L + SoilWeight(L, D, t, gamma))/FoSsoil                   
        Vmax = min(Vmax1, Vmax2, Vmax3)
               
    elif soil_type == 'sand':
        # Setting default gamma values per soil type [kN/m3]
        gamma = 9
        # Definition for non-cohesive soils
        Nq = np.e**(np.pi*np.tan(np.radians(phi)))*np.tan(np.radians(45) + np.radians(phi)/2)**2 # Lateral-bearing capacity factor
        sigma_av_L = gamma*2*L/3                             # Effective stress (average)
        Hmax = 0.5*D*Nq*gamma*L**2

        M = - resultsLoad['V']*rlugTilt(rlug,zlug,thetalug) - resultsLoad['H']*(zlugTilt(rlug,zlug,thetalug) - zlug)
        
        # Torsion capacity
        Fo = PileSurface(L, D)*beta*sigma_av_L
        To = Fo; print(To)
        Ti = Fo
        Tmax = Ti + To
        
        # Introduce twist effects due to installation misaligment
        T = resultsLoad['H']*rlug*np.sin(np.deg2rad(psilug))
        nhuT = T/Tmax; nhuV = resultsLoad['H']/Fo;
        nhuVstar = np.sqrt(nhuV**2 - nhuT**2)
        alphastar = alpha*(nhuVstar/nhuV)
    
        # "Coring"        
        Vmax = PileWeight(L, D, t, rhows) + PileSurface(L, D)*beta*sigma_av_L + PileSurface(L,(D - 2*t))*beta*sigma_av_L
        def y(depth):
            return np.e**(-depth) - 1 + depth
        Ze = D/(4*7); Zi = D/(4*5)    
        Vmax = 7*gamma*Ze**2*y(L/Ze)*PileSurface(L, D)/L + 5*gamma*Zi**2*y(L/Zi)*PileSurface(L,(D - 2*t))/L
        # Vmax = PileWeight(L, D, t, rhows) + gamma*L**2/(2*(beta + beta)*np.pi*D)
     
    # Submerged pile weight (inc. stiffening plus vent) assessed as a factor
    Wp = 1.00*PileWeight(L, D, t, (rhows)) 
    # Submerged weight of the soil plug
    Ws = SoilWeight(L, D, t, gamma)
    
    # Capacity envelope
    aVH = 0.5 + lambdap; bVH = 4.5 + lambdap/3 
    print('Env. exp = ' +str(aVH)+'   '+str(bVH))
    UC = (resultsLoad['H']/Hmax)**aVH + (resultsLoad['V']/Vmax)**bVH      
    x = np.cos(np.linspace (0, np.pi/2, 100))
    y = (1 - x**bVH)**(1/aVH)
    X = Hmax*x; Y = Vmax*y
    plt.plot(X, Y, color = 'b')
    plt.scatter(resultsLoad['H'], resultsLoad['V'], color = 'r')
        
    # Set labels and title
    plt.xlabel('Horizontal capacity [kN]')
    plt.ylabel('Vertical capacity [kN]')
    plt.suptitle('VH suction pile capacity envelope')
    plt.axis([0,1.3*max(X[0],resultsLoad['H']),0,1.3*max(Y[-1],resultsLoad['V'])]) 
    plt.grid(True)
    plt.show()
    
    resultsSuction = {}
    if soil_type == 'clay':
        resultsSuction['Horizontal max.'] = Hmax[0] # Capacity at specified loading angle in clay
    elif soil_type == 'sand':   
        resultsSuction['Horizontal max.'] = Hmax    # Capacity at specified loading angle in sand
    resultsSuction['Vertical max.'] = Vmax          # Capacity at specified loading angle
    if soil_type == 'clay':
        resultsSuction['UC'] = UC[0]                # Unity check in clay
    elif soil_type == 'sand':
        resultsSuction['UC'] = UC                   # Unity check in sand
    resultsSuction['Weight Pile'] = Wp              # in kN
    resultsSuction['Weight Soil'] = Ws              # in kN
    resultsSuction['t'] = t                         # in m
    
    return resultsSuction
        
def getCapacitySuctionSimp(D, L, Su0=15.0, k=2.0, alpha=0.515, rhows=66.90):

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
    safety_factor: float
        Specify 'yes' or 'no'. This affects load and resistance factors.    
    Su0 : float 
        Undrained shear strength at the mudline (clay only), default is 2.39. [kPa]
    k : float 
        Undrained shear strength gradient (clay only), default is 1.41. [kPa/m]
    alpha : float 
        Skin friction coefficient (outer and inner - clay only), default is 0.7 [-] 
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
   
    gamma = 8.75                            # Effective unit weight of clay
    lambdap = L/D;                       # Suction pile slenderness ratio   
    t = 10*D/1e3                         # Thickness of the pile
    Np_fixed = 10.25;                    # From Np vs L/D chart from CAISSON_VHM
    Su_av_L = Su0 + k*zlug;              # Undrained shear strength values (average)
    print(Su_av_L)
    Su_tip = Su0 + k*L;                  # Undrained shear strength values (tip)
    print(Su_tip)
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
    
    Hmax = Np_fixed*L*D*Su_av_L/FoSh
    # "Plugged" (Reverse end bearing capacity - passive suction) 
    Vmax1 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alpha*Su_av_L + Nc*Su_tip*(np.pi/4)*D**2)/FoSsoil
    # "Coring"        
    Vmax2 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alpha*Su_av_L + PileSurface(L,(D - 2*t))*alpha*Su_av_L)/FoSsoil
    # "Leaking"        
    Vmax3 = (PileWeight(L, D, t, rhows) + PileSurface(L, D)*alpha*Su_av_L + SoilWeight(L, D, t, gamma))/FoSsoil 
    print(Vmax1); print(Vmax2); print(Vmax3)                  
    Vmax = min(Vmax1,Vmax2,Vmax3)
    
    # Submerged pile weight (inc. stiffening plus vent) assessed as a factor
    Wp = 1.00*PileWeight(L, D, t, (rhows)) 
    # Submerged weight of the soil plug
    Ws = SoilWeight(L, D, t, gamma) 
    
    H = Tm*np.cos(np.deg2rad(thetam)); V = Tm*np.sin(np.deg2rad(thetam))
    
    # Capacity envelope
    aVH = 0.5 + lambdap; bVH = 4.5 + lambdap/3 
    print('Env. exp =' +str(aVH)+str(bVH))
    UC = (H/Hmax)**aVH + (V/Vmax)**bVH

    x = np.cos(np.linspace (0,np.pi/2,1000))
    y = (1 - x**bVH)**(1/aVH)
    X = Hmax*x; Y = Vmax*y

    iload = np.argwhere(np.diff(np.sign(np.degrees(np.arctan(Y/X)) - thetam))).flatten()[0]
    #iload = np.argwhere(np.diff(np.sign(np.degrees(np.arctan(Y/X)) - thetam_sf))).flatten()[0]
    H_good = X[iload]
    V_good = Y[iload]

    #H_good = Hmax*np.exp(np.log(0.5)/aVH)
    #V_good = Vmax*np.exp(np.log(0.5)/bVH)

    if plot:
        plt.plot(X,Y,color = 'c')
        plt.scatter(H,V,color = 'y')
        plt.scatter(H_good, V_good, color='g')
        # Set labels and title
        plt.xlabel('Horizontal capacity [kN]')
        plt.ylabel('Vertical capacity [kN]')
        plt.suptitle('VH suction pile capacity envelope SIMP')
        plt.axis([0,1.3*max(X[0],H),0,1.3*max(Y[-1],V)]) 
        plt.grid(True)
        plt.show()
    
    resultsSuctionSimp = {}
    resultsSuctionSimp['Horizontal max.'] = Hmax    # Capacity at specified loading angle
    resultsSuctionSimp['Vertical max.'] = Vmax      # Capacity at specified loading angle
    resultsSuctionSimp['UC'] = UC                   # Unity check
    resultsSuctionSimp['Weight Pile'] = Wp          # in kN
    resultsSuctionSimp['Weight Soil'] = Ws          # in kN
    resultsSuctionSimp['t'] = t
    resultsSuctionSimp['H_good'] = H_good
    resultsSuctionSimp['V_good'] = V_good
    
    return resultsSuctionSimp



if __name__ == '__main__':
          

    D = 5.45
    L = 16.4
    zlug = 9.32

    Tm = 23000
    thetam = 5
    

    resultsSuction = getCapacitySuction(D, L, Tm, thetam, zlug, safety_factor='no')

    resultsSuctionSimp = getCapacitySuctionSimp(D, L, Tm, thetam, zlug, safety_factor='no')
    
    print('*************** Suction Pile Result Simp *********************')

    print('Anchor thickness,                    ' , resultsSuctionSimp['t'], '[m]')
    print('Anchor steel weight,                 ' , resultsSuctionSimp['Weight Pile'], '[kN]')
    print('Soil plug weight,                    ' , resultsSuctionSimp['Weight Soil'], '[kN]')
    print('Horizontal max. capacity,            ' , resultsSuctionSimp['Horizontal max.'], '[kN]')
    print('Vertical max. capacity,              ' , resultsSuctionSimp['Vertical max.'], '[kN]') 
    print('Unity check capacity,                ' , resultsSuctionSimp['UC'], '[-]') 

    print('**************************************************************') 



    ''' 
    Testing the function 
    
    # Retrieves input data from a separate config file
    with open('SuctionConfig.yml', 'r') as f:
        configSuction = yaml.full_load(f)
    D = configSuction['D']; L = configSuction['L'];  
    thetalug = configSuction['thetalug']; psilug = configSuction['psilug'];
    soil_type = configSuction['soil_type'];
    Su0 = configSuction['Su0']; k = configSuction['k'];    
    alpha = configSuction['alpha']; beta = configSuction['beta'];
    rhows = configSuction['rhows'];
   
    # Retrieves input data from a separate config file
    with open('LoadConfig_FMO.yml', 'r') as f:
        configLoad_FMO = yaml.full_load(f)
    Tm = configLoad_FMO ['Tm']; thetam = configLoad_FMO ['thetam']
    zlug = configLoad_FMO ['zlug']; 
    soil_type = configLoad_FMO ['soil_type'];
    line_type = configLoad_FMO ['line_type'];
    Su0 = configLoad_FMO ['Su0']; k = configLoad_FMO['k']
    d = configLoad_FMO ['d']; 
   
   '''
    
    print('******************  Suction Pile Result  *********************')

    print('Anchor thickness,                    ' , resultsSuction['t'], '[m]')
    print('Anchor steel weight,                 ' , resultsSuction['Weight Pile'], '[kN]') 
    print('Soil plug weight,                    ' , resultsSuction['Weight Soil'], '[kN]')
    print('Horizontal max. capacity,            ' , resultsSuction['Horizontal max.'], '[kN]')
    print('Vertical max. capacity,              ' , resultsSuction['Vertical max.'], '[kN]') 
    print('Unity check capacity,                ' , resultsSuction['UC'], '[-]') 

    print('**************************************************************') 
          
    