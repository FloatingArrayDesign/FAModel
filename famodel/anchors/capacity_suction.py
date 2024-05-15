# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:53:52 2024

@author: fmoreno
"""

import yaml      # Allow access to config file for user inputs
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from capacity_load import getAnchorLoad

def getCapacitySuction(D,L,thetalug,psilug,Su0,dSu,gamma, 
                    alphao,alphai,Np,nhu,Ab,rhows):
    '''Calculate the inclined load capacity of a Suction Caisson Anchor in sand or clay.
    The calculation is based on the soil properties, anchor geometry, and the angle of inclined load.  
    
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
    Su0 : float 
        Undrained shear strength at the mudline (clay only), default is 2.39. [kPa]
    dSu : float 
        Undrained shear strength gradient (clay only), default is 1.41. [kPa/m]
    gamma: float 
        Effective unit weight of the sand soil, default is 9 for sand, 4.7 for clay. [kN/m3] 
    alphao : float 
        Skin friction coefficient (outer), default is 0.7 [-] 
    alphai : float
        Skin friction coefficient (inner), default is 0.7 [-]
    Np : float
        Friction angle of the sand soil (sand only), default is 30. [deg]  
    nhu : float 
        The API- Factor [---], default is 0.5. (sand only)  
    rhows : float
        Submerged steel density [t/m3]
    
    Returns
    -------
    Hmax : float 
        Maximum horizontal capacity [kN]
    Vmax : float 
        Maximum vertical capacity [kN]
    '''
    
    lambdap = L/D; m = -2/3;
    t = 10*D/1e3          # Thickness of the pile
    Nc = min (6.2*(1 + 0.34*np.arctan(lambdap)),9)
    ez = (Su0*L**2/2 + dSu*L**3/3)/(Su0*L + dSu*L**2/2)
    Np_fixed = 10.5; Np_free = 4 # From Np vs L/D chart from CAISSON_VHM
    Su_av_L = Su0 + dSu*L/3; Su_tip = Su0 + dSu*L
    zlug = ez           # Optimized depth of the lug 
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
    
    resultsLoad = getAnchorLoad(Tm,thetam,zlug,Su0,dSu,nhu,Ab,Nc)
    M = - resultsLoad['V']*rlugTilt(rlug,zlug,thetalug) - resultsLoad['H']*(zlugTilt(rlug,zlug,thetalug) - ez)
    
    # M modifies the Hmax capacity
    def f(Hmax):
         return m*(Hmax/(L*D*(Su0 + dSu*L/3)) - Np_fixed) + M*(Hmax/(L*D*(Su0 + dSu*L/3))/(Hmax*L))
    Hmax = fsolve(f,5); print(Hmax[0])
    
    Fo = PileSurface(L,D)*alphao*Su_av_L
    T = resultsLoad['H']*rlug*np.sin(np.deg2rad(psilug))
    To = PileSurface(L,D)*alphao*Su_av_L
    Ti = PileSurface(L,(D -2*t))*alphai*Su_av_L
    Tbase = np.pi*D**3*Su_tip/12
    Tmax = min(Ti + To,To + Tbase)
    
    # Introduce twist effects
    nhuT = T/Tmax; nhuV = resultsLoad['H']/Fo;
    nhuVstar = np.sqrt(nhuV**2 - nhuT**2)
    alphaostar = alphao*(nhuVstar/nhuV)
    
    Vmax1 = PileWeight(L,D,t,rhows) + PileSurface(L,D)*alphaostar*Su_av_L + Nc*Su_tip*np.pi*D**2                    # "Plugged" (Reverse end bearing capacity) 
    Vmax2 = PileWeight(L,D,t,rhows) + PileSurface(L,D)*alphaostar*Su_av_L + PileSurface(L,(D - 2*t))*alphai*Su_av_L # "Coring"
    Vmax3 = PileWeight(L,D,t,rhows) + PileSurface(L,D)*alphaostar*Su_av_L + SoilWeight(L,D,gamma)                   # "Leaking"
    Vmax = min(Vmax1,Vmax2,Vmax3)
    
    print('Hmax = ' +str(Hmax[0]) +' kN')
    print('Vmax = ' +str(Vmax) +' kN')
    print('Tmax = ' +str(Tmax) +' kN*m')
    
    Wp = 1.10*PileWeight(L,D,t,rhows)  
    
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
    
if __name__ == '__main__':
          
    ''' 
    Testing the function 
    '''
  
    # Retrieves input data from a separate config file
    with open('SuctionConfig.yml', 'r') as f:
        configSuction = yaml.full_load(f)
    D = configSuction['D']; L = configSuction['L'];  
    thetalug = configSuction['thetalug']; psilug = configSuction['psilug'];
    Su0 = configSuction['Su0']; dSu = configSuction['dSu']; gamma = configSuction['gamma']
    alphao = configSuction['alphao']; alphai = configSuction['alphai']
    Np = configSuction['Np']; Nc = configSuction['Nc']; Ab = configSuction['Ab']
    nhu = configSuction['nhu']; rhows = configSuction['rhows']
    
    # Retrieves input data from a separate config file
    with open('LoadConfig.yml', 'r') as f:
        configLoad = yaml.full_load(f)
    Tm = configLoad['Tm']; thetam = configLoad['thetam']
    zlug = configLoad['zlug']; 
    Su0 = configLoad['Su0']; dSu = configLoad['dSu']
    nhu = configLoad['nhu']; 
    Nc = configLoad['Nc']; Ab = configLoad['Ab']  
    
    resultsSuction = getCapacitySuction(D,L,thetalug,psilug,Su0,dSu,gamma, 
                        alphao,alphai,Np,nhu,Ab,rhows)
    
    print('******************  Suction Pile Result  *********************')

    print('Anchor thickness,                    ' , resultsSuction['t'], '[m]')
    print('Anchor steel weight,                 ' , resultsSuction['Weight'], '[t]') 
    print('Horizontal max. capacity,            ' , resultsSuction['Horizontal max.'], '[kN]')
    print('Vertical max. capacity,              ' , resultsSuction['Vertical max.'], '[kN]') 
    print('Unity check capacity,                ' , resultsSuction['UC'], '[-]') 

    print('**************************************************************') 
     
    # delta_phiMH = 2.5 - 0.25*lambdap  # Formula in deg, depends on ez_L
    # phiMH = -np.arctan(ez_L) - np.deg2rad(delta_phiMH)
    # aMH = np.cos(phiMH)/Np_fixed; bMH = Np_free*np.sin(phiMH)
    # print(aMH); print(bMH)
    
    # N = 50
    # alpha = np.linspace(0,2*np.pi,N)
    # phi = np.linspace(0,np.pi/2,N)
    # Alpha, Phi = np.meshgrid(alpha,phi)
    
    # x = (aMH*np.cos(phiMH)*np.cos(Alpha) - bMH*np.sin(phiMH)*np.sin(Alpha))*np.sin(Phi)
    # y = (aMH*np.sin(phiMH)*np.cos(Alpha) + bMH*np.cos(phiMH)*np.sin(Alpha))*np.sin(Phi)
    # z = ((1-np.abs(1/Vmax))**bVH)**(2/aVH)*np.cos(Phi)
    
    # # First subplot
    # fig = plt.figure(figsize=(15,20))
    # ax = fig.add_subplot(3,1,1,projection='3d')
    # surf = ax.plot_surface(x,y,z,cmap=cm.coolwarm)
    # load = ax.scatter(H/(D*L*(Su0+dSu*L/3)), M/(D*L*L*(Su0+dSu*L/3)), V/Vmax, color = 'r')
    
    # # Set labels and title
    # ax.set_xlabel('H* (H/(D*L*Su,av,L)) [-]')
    # ax.set_ylabel('M* (M/(D*L^2*Su,av,L)) [-]')
    # ax.set_zlabel('V [-]')
    # ax.set_title('VHM Ellipsoid suction pile capacity')
    # plt.show()