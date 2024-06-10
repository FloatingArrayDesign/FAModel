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


def getCapacitySuction(D, L, Tm, thetam, zlug, safety_factor='yes', line_type='chain', d='0.15', thetalug =5, psilug=7.5, soil_type='clay', 
                Su0=10.0, k=2.0, alpha=0.7, beta=0.46, rhows=6.85):
    '''Calculate the inclined load capacity of a suction pile in sand or clay.
    The calculation is based on the soil properties, anchor geometry and inclineded load.  
    
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
        Undrained shear strength at the mudline (clay only), default is 2.39. [kPa]
    k : float 
        Undrained shear strength gradient (clay only), default is 1.41. [kPa/m]
    alpha : float 
        Skin friction coefficient (outer and inner - clay only), default is 0.7 [-] 
    beta : float
        Skin friction coefficient (sand only), default is 0.46 [-]
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
    
    # Setting default safety factors
    # API RP 2GEO reverse end bearing capacity mimumum recommended. [-]
    if safety_factor == 'yes':
        FoSsoil1 = 2.0
        FoSsoil2 = 1.5
    elif safety_factor == 'no':
        FoSsoil1 = 1.0
        FoSsoil2 = 1.0
          
    # Setting default gamma values per soil type
    # Effective unit weight of the soil, default 17 for clay. [kN/m3]
    if soil_type == 'clay':
        gamma = 17
    elif soil_type == 'sand':
        gamma = 9
    
    lambdap = L/D; m = -2/3;                        # Suction pile slenderness ratio
    t = 5*D/1e3                                     # Thickness of the pile
    Nc = min (6.2*(1 + 0.34*np.arctan(lambdap)),9)  # End-bearing capacity factor
    ez = (Su0*L**2/2 + k*L**3/3)/(Su0*L + k*L**2/2)
    Np_fixed = 10.5; Np_free = 4                    # From Np vs L/D chart from CAISSON_VHM
    Su_av_L = Su0 + k*L/3; Su_tip = Su0 + k*L       # Undrained shear strength values (average, tip)
    sigma_av_L = gamma*2*L/3                        # Effective stress (average, tip)
    #zlug = ez          # Optimized depth of the lug 
    rlug = D/2          # Radial position of the lug
    
    # Outer and inner surface of the pile skirt
    def PileSurface(Len,Dia):
        Sp = np.pi*Dia*Len
        return Sp    
    # Dry and wet mass of the pile    
    def PileWeight(Len,Dia,tw,rho):
        Wp = ((np.pi/4)*(Dia**2-(Dia-tw)**2)*Len + np.pi*Dia**2*tw)*rho
        return Wp 
    # Mass of the soil plug      
    def SoilWeight(Len,Dia,gamma_soil): 
        Wsoil =(np.pi/4)*Dia**2*Len*gamma_soil
        return Wsoil
    # Tilt and twist effects due to installation misaligments
    def rlugTilt(r,z,theta):
        R = r*np.cos(np.deg2rad(theta)) - z*np.sin(np.deg2rad(theta))
        return R   
    def zlugTilt(r,z,theta):
        Z = r*np.sin(np.deg2rad(theta)) + z*np.cos(np.deg2rad(theta))
        return Z
           
    Hmax = Np_fixed*L*D*Su_av_L; H0 = Np_free*L*D*Su_av_L;
    Mmax = Np_fixed*L*L*D*Su_av_L; print(Hmax)
    
    # Introduce tilt effects due to installation misaligment
    resultsLoad = getAnchorLoadDNV(Tm=5000, thetam=25, zlug=14)
    M = - resultsLoad['V']*rlugTilt(rlug,zlug,thetalug) - resultsLoad['H']*(zlugTilt(rlug,zlug,thetalug) - ez)
    
    # M modifies the Hmax capacity
    def f(Hmax):
         return m*(Hmax/(L*D*(Su0 + k*L/3)) - Np_fixed) + M*(Hmax/(L*D*(Su0 + k*L/3))/(Hmax*L))
    Hmax = fsolve(f,5); print(Hmax[0])
    
    # Torsion capacity
    if soil_type == 'clay':
        Fo = PileSurface(L,D)*alpha*Su_av_L
        To = Fo
        Ti = PileSurface(L,(D -2*t))*alpha*Su_av_L
        Tbase = np.pi*D**3*Su_tip/12
        Tmax = min(Ti + To,To + Tbase)        
    elif soil_type == 'sand':
        Fo = PileSurface(L,D)*beta*sigma_av_L
        To = Fo
        Ti = Fo
        Tmax = Ti + To      
    
    # Introduce twist effects due to installation misaligment
    T = resultsLoad['H']*rlug*np.sin(np.deg2rad(psilug))
    nhuT = T/Tmax; nhuV = resultsLoad['H']/Fo;
    nhuVstar = np.sqrt(nhuV**2 - nhuT**2)
    alphastar = alpha*(nhuVstar/nhuV)
    
    if soil_type == 'clay':
        # "Plugged" (Reverse end bearing capacity - passive suction) 
        Vmax1 = (PileWeight(L,D,t,rhows) + PileSurface(L,D)*alphastar*Su_av_L + Nc*Su_tip*np.pi*D**2)/FoSsoil1
        # "Coring"        
        Vmax2 = (PileWeight(L,D,t,rhows) + PileSurface(L,D)*alphastar*Su_av_L + PileSurface(L,(D - 2*t))*alpha*Su_av_L)/FoSsoil2
        # "Leaking"        
        Vmax3 = (PileWeight(L,D,t,rhows) + PileSurface(L,D)*alphastar*Su_av_L + SoilWeight(L,D,gamma))/FoSsoil2                   
        Vmax = max(Vmax1,Vmax2,Vmax3)
    elif soil_type == 'sand':
        # "Coring"        
        Vmax = PileWeight(L,D,t,rhows) + PileSurface(L,D)*beta*sigma_av_L + PileSurface(L,(D - 2*t))*beta*sigma_av_L
    
    print('Hmax = ' +str(Hmax[0]) +' kN')
    print('Vmax = ' +str(Vmax) +' kN')
    print('Tmax = ' +str(Tmax) +' kN*m')
    
    # Top plate weight (inc. stiffening plus vent) assessed as a factor
    Wp = 1.15*PileWeight(L,D,t,(rhows + 1))  
    
    # Capacity envelope
    aVH = 0.5 + lambdap; bVH = 4.5 + lambdap/3 
    UC = (resultsLoad['H']/Hmax)**aVH + (resultsLoad['V']/Vmax)**bVH      
    x = np.cos(np.linspace (0,np.pi/2,100))
    y = (1 - x**bVH)**(1/aVH)
    X = Hmax*x; Y = Vmax*y
    plt.plot(X,Y,color = 'b')
    plt.scatter(resultsLoad['H'],resultsLoad['V'],color = 'r')
    
    # Set labels and title
    plt.xlabel('Horizontal capacity [kN]')
    plt.ylabel('Vertical capacity [kN]')
    plt.suptitle('VH suction pile capacity envelope')
    plt.axis([0,1.3*max(X[0],resultsLoad['H']),0,1.3*max(Y[-1],resultsLoad['V'])]) 
    plt.grid(True)
    plt.show()
    
    resultsSuction = {}
    resultsSuction['Horizontal max.'] = Hmax[0] # Capacity at specified loading angle
    resultsSuction['Vertical max.'] = Vmax      # Capacity at specified loading angle
    resultsSuction['UC'] = UC[0]                # Unity check
    resultsSuction['Weight'] = Wp
    resultsSuction['t'] = t
    
    return resultsSuction
    


     
def getCapacitySuctionSimp(D, L, Tm, thetam, zlug, safety_factor='yes', Su0=10.0, k=2.0, alpha=0.7, rhows=6.85, plot=True):
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
    
    # Setting default safety factors
    # API RP 2SK permanent intact condition lateral/axial safety factors. [-]
    # API RP 2GEO reverse end bearing capacity mimumum recommended. [-]
    if safety_factor == 'yes':
        FoSh = 1.6
        FoSv = 2.0
        FoSsoil = 2.0
    elif safety_factor == 'no':
        FoSh = 1.0
        FoSv = 1.0
        FoSsoil = 1.0
    
    gamma = 4.7                          # Effective unit weight of clay
    lambdap = L/D;                       # Suction pile slenderness ratio   
    t = 5*D/1e3                          # Thickness of the pile
    Np_fixed = 10.5;                     # From Np vs L/D chart from CAISSON_VHM
    Su_av_L = Su0 + k*L/3;               # Undrained shear strength values (average)
    Su_tip = Su0 + k*L;                  # Undrained shear strength values (tip)
    Nc = min (6.2*(1 + 0.34*np.arctan(lambdap)),9)  # End-bearing capacity factor
    
    # Outer and inner surface of the pile skirt
    def PileSurface(Len,Dia):
        Sp = np.pi*Dia*Len
        return Sp    
    # Dry and wet mass of the pile    
    def PileWeight(Len,Dia,tw,rho):
        Wp = ((np.pi/4)*(Dia**2-(Dia-tw)**2)*Len + np.pi*Dia**2*tw)*rho
        return Wp 
    # Mass of the soil plug      
    def SoilWeight(Len,Dia,gamma_soil): 
        Wsoil =(np.pi/4)*Dia**2*Len*gamma_soil
        return Wsoil
    
    Hmax = Np_fixed*L*D*Su_av_L / FoSh
    Vmax = (PileWeight(L,D,t,rhows) + PileSurface(L,D)*alpha*Su_av_L + Nc*Su_tip*np.pi*D**2) / FoSsoil / FoSv
    
    Wp = 1.15*PileWeight(L,D,t,(rhows + 1)) 
    
    H = Tm*np.cos(np.deg2rad(thetam)); V = Tm*np.sin(np.deg2rad(thetam))
    #H = FoSh*Tm*np.cos(np.deg2rad(thetam)); V = FoSv*Tm*np.sin(np.deg2rad(thetam))
    #thetam_sf = np.degrees(np.arctan(V_sf/H_sf))
    
    # Capacity envelope
    aVH = 0.5 + lambdap; bVH = 4.5 + lambdap/3 
    UC = (H/Hmax)**aVH + (V/Vmax)**bVH
    #UC_sf = (H_sf/Hmax)**aVH + (V_sf/Vmax)**bVH

    #print(Hmax, Vmax)
    
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
        plt.plot(X,Y,color = 'b')
        plt.scatter(H,V,color = 'r')
        #plt.scatter(H_sf,V_sf,color = 'm')
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
    resultsSuctionSimp['Weight'] = Wp
    resultsSuctionSimp['t'] = t
    resultsSuctionSimp['H_good'] = H_good
    resultsSuctionSimp['V_good'] = V_good
    
    return resultsSuctionSimp



if __name__ == '__main__':
          

    D = 10
    L = 25
    D = 3.02656557
    L = 8.95993379

    #fx = 54735571.30174828
    #fy = 42519484.98012456
    fx = 4522222.788895202
    fy = 2948278.926831712

    Tm = np.linalg.norm([fx, fy])/1000
    thetam = np.degrees(np.arctan(fy/fx))
    zlug = 2.0

    #resultsSuction = getCapacitySuction(D, L, Tm, thetam, zlug, safety_factor='no')

    resultsSuctionSimp = getCapacitySuctionSimp(D, L, Tm, thetam, zlug, safety_factor='yes')
    
    print('*************** Suction Pile Result Simp *********************')

    print('Anchor thickness,                    ' , resultsSuctionSimp['t'], '[m]')
    print('Anchor steel weight,                 ' , resultsSuctionSimp['Weight'], '[t]') 
    print('Horizontal max. capacity,            ' , resultsSuctionSimp['Horizontal max.'], '[kN]')
    print('Vertical max. capacity,              ' , resultsSuctionSimp['Vertical max.'], '[kN]') 
    print('Unity check capacity,                ' , resultsSuctionSimp['UC'], '[-]') 

    print('**************************************************************') 



    ''' 
    Testing the function 
    '''
    """
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
    
    resultsSuction = getCapacitySuction(D, L, Tm, thetam, zlug, line_type, d, thetalug =5, psilug=7.5, soil_type='clay', 
                    Su0=2.39, k=1.41, alpha=0.7, beta=0.46, rhows=6.85)
    
    print('******************  Suction Pile Result  *********************')

    print('Anchor thickness,                    ' , resultsSuction['t'], '[m]')
    print('Anchor steel weight,                 ' , resultsSuction['Weight'], '[t]') 
    print('Horizontal max. capacity,            ' , resultsSuction['Horizontal max.'], '[kN]')
    print('Vertical max. capacity,              ' , resultsSuction['Vertical max.'], '[kN]') 
    print('Unity check capacity,                ' , resultsSuction['UC'], '[-]') 

    print('**************************************************************') 
          
    """