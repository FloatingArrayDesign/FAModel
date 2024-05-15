"""Torpedo anchor capacity calculation functions in 'level 2' model, 
calculating inclined load capacity over the range of angles.
Lead author: Felipe Moreno. 
"""

import yaml      # Allow access to config file for user inputs
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from matplotlib import cm
from mpl_toolkits import mplot3d 
from capacity_load import getAnchorLoad
   
def getCapacityTorpedo(D1,D2,L1,L2,Hp,Su0,dSu,gamma,alphao,Np,rhows):
    
    '''Calculate the inclined load capacity of a Suction Caisson Anchor in sand or clay.
    The calculation is based on the soil properties, anchor geometry, and the angle of inclined load.  
    
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
    Hp : float 
        Torpedo pile embeded depth. [m]        
    A_angle : float 
        The angle of the tension force at the padeye, default is 60. [deg] 
    F_angle : float 
        The adjustment factor of the load angle at the padeye, default is 2. [deg] 
    soil_type : string
        Specify 'sand' or 'clay'. This affects what other soil parameters are used.
    phi : float
        The friction angle of the sand soil (sand only), default is 30. [deg] 
    gamma: float 
        The effective unit weight of the sand soil, default is 9 for sand, 4.7 for clay. [kN/m3] 
    J: float 
        The API- Factor [---], default is 0.5. (sand only)
    Su0 : float 
        The Undrained shear strength at the mudline (clay only), default is 2.39. [kPa]
    k : float 
        The Undrained shear strength gradient (clay only), default is 1.41. [kPa/m]
    alpha : float
        The skin friction coefficient (clay only), default is 0.7. 
    rhows : float
        Submerged steel density [t/m3]        
    
    Returns
    -------
    Ft_values : flaot 
        The inlcind force capacity at the padeye [kN]
    Anchor_vol : floar 
        The corresponding anchor steel volume [m3]
    '''

    Hp = 10; m = -2/3;
    L = L1  + L2;
    D = (D1*L1 +D2*L2)/L;     
    rlug = D2/2; zlug = 0.75*L; 
    t = 100*D/1e3
    a = Hp; b = Hp + L1; c = Hp + L1 + L2;
    Ballast = 10     # tonnes of ballast to lower the CoG
    
    ez_Su_den = D1*Su0*(b - a) + 0.5*D1*dSu*(b**2 - a**2) + D2*Su0*(c - b) - 0.5*D2*dSu*(c**2 - b**2)
    ez_Su_num = D1*Su0*(a**2 - a*b) + 0.33*D1*dSu*(b**3 - a**3) + b**2*(0.5*D1*Su0 - 0.5*D1*a*dSu) - a**2*(0.5*D1*Su0 - 0.5*D1*a*dSu)\
        + D2*Su0*(b**2 - b*c) + 0.33*D2*dSu*(c**3 - b**3) + c**2*(0.5*D2*Su0 - 0.5*D2*b*dSu) - b**2*(0.5*D2*Su0 - 0.5*D2*b*dSu)
    ez_Su = ez_Su_num/ez_Su_den
    ez_Su_L = ez_Su/L
    Np_fixed = 10; Np_free = 4 # From Np vs L/D chart from CAISSON_VHM
       
    # Dry and wet mass of the pile
    def PileSurface(Len1,Len2,Dia1,Dia2):
        Sp = np.pi*Dia1*Len1 + 4*Len2*(Dia2 - Dia1)
        return(Sp)  
    # Dry and wet mass of the pile    
    def PileWeight(Len1,Len2,Dia1,Dia2,tw,rho):
        Wp = ((np.pi/4)*(Dia1**2-(Dia1-tw)**2)*(Len1+Len2) + 4*Len2*Dia2*tw)*rho
        return(Wp)
           
    Hmax = L*D*Np*(Su0 + dSu*Hp + dSu*L/3); Mmax = Hmax*(ez_Su_L - zlug)
    print('Hmax = ' +str(Hmax))
    Vmax = PileWeight(L1,L2,D1,D2,t,rhows) + PileSurface(L1,L2,D1,D2)*alphao*(Su0 + dSu*Hp + dSu*L/3) + Ballast*9.81
    print('Vmax = ' +str(Vmax))
    
    resultsLoad = getAnchorLoad(Tm,thetam,zlug,Su0,dSu,nhu,Ab,Nc)
    M = - resultsLoad['V']*rlug - resultsLoad['H']*(zlug - ez_Su)
    
    # M modifies the Hmax capacity
    def f(Hmax):
         return m*(Hmax/(L*D*(Su0 + dSu*Hp + dSu*L/3)) - Np_fixed) + M*(Hmax/(L*D*(Su0 + dSu*Hp + dSu*L/3))/(Hmax*L))
    Hmax = fsolve(f,5); print(Hmax[0])
    
    Ta = Tm/(np.e**(nhu*(resultsLoad['angle'] - thetam)))
    aVH = 0.5 + L/D; bVH = 4.5 + L/(3*D) 
    H = Ta*np.cos(np.deg2rad(resultsLoad['angle'])); V = Ta*np.sin(np.deg2rad(resultsLoad['angle']))
    UC = (H/Hmax)**aVH + (V/Vmax)**bVH
       
    x = np.cos(np.linspace (0,np.pi/2,100))
    y = (1 - x**bVH)**(1/aVH)
    X = Hmax*x; Y = Vmax*y
    plt.plot(X,Y,color = 'b')
    plt.scatter(H,V,color = 'r')
       
    # Set labels and title
    plt.xlabel('Horizontal capacity [kN]')
    plt.ylabel('Vertical capacity [kN]')
    plt.suptitle('VH torpedo pile capacity envelope')
    plt.axis([0,1.3*max(X[0],resultsLoad['H']),0,1.3*max(Y[-1],resultsLoad['V'])]) 
    plt.grid(True)
    plt.show()
     
    delta_phiMH = 2.5 - 0.25*(L/D)  # Formula in deg, depends on ez_su_L
    phiMH = -np.arctan(ez_Su_L) - np.deg2rad(delta_phiMH)
    aMH = Np_fixed/np.cos(phiMH); bMH = -Np_free*np.sin(phiMH)
       
    N = 50
    alpha = np.linspace(0, 2*np.pi, N)
    phi = np.linspace(0, np.pi/2, N)
    Alpha, Phi = np.meshgrid(alpha, phi)
       
    x = (aMH*np.cos(phiMH)*np.cos(Alpha) - bMH*np.sin(phiMH)*np.sin(Alpha))*np.sin(Phi)
    y = (aMH*np.sin(phiMH)*np.cos(Alpha) + bMH*np.cos(phiMH)*np.sin(Alpha))*np.sin(Phi)
    z = ((1-np.abs(1/Vmax))**bVH)**(2/aVH)*np.cos(Phi)
       
    # First subplot
    fig = plt.figure(figsize=(15,20))
    ax = fig.add_subplot(3,1,1, projection='3d'); #ax.view_init(10,90)
    surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm)
    # load = ax.scatter(Hmax/(D*L*(Su0+dSu*L/3)), Mmax/(D*L*L*(Su0+dSu*L/3)), V/Vmax, color = 'r')
       
    # Set labels and title
    ax.set_xlabel('H* (H/(D*L*Su,av,L)) [-]')
    ax.set_ylabel('M* (M/(D*L^2*Su,av,L)) [-]')
    ax.set_zlabel('V [-]')
    ax.set_title('VHM Ellipsoid torpedo pile capacity')
    plt.show() 
    
    resultsTorpedo = {}
    resultsTorpedo['Horizontal max.'] = Hmax[0]    # Capacity at specified loading angle
    resultsTorpedo['Vertical max.'] = Vmax      # Capacity at specified loading angle
    resultsTorpedo['UC'] = UC[0]                   # Unity check
   
    return resultsTorpedo
    
if __name__ == '__main__':
          
    ''' 
    Testing the function 
    '''
  
    # Retrieves input data from a separate config file
    with open('TorpedoConfig.yml', 'r') as f:
        configTorpedo = yaml.full_load(f)
    D1 = configTorpedo['D1']; D2 = configTorpedo['D2'];
    L1 = configTorpedo['L1']; L2 = configTorpedo['L2']; 
    Hp = configTorpedo['Hp']
    Su0 = configTorpedo['Su0']; dSu = configTorpedo['dSu']; gamma = configTorpedo['gamma']
    alphao = configTorpedo['alphao']; 
    Np = configTorpedo['Np']; rhows = configTorpedo['rhows'] 
    
    # Retrieves input data from a separate config file
    with open('LoadConfig.yml', 'r') as f:
        configLoad = yaml.full_load(f)
    Tm = configLoad['Tm']; thetam = configLoad['thetam']
    zlug = configLoad['zlug']; 
    Su0 = configLoad['Su0']; dSu = configLoad['dSu']
    nhu = configLoad['nhu']; 
    Nc = configLoad['Nc']; Ab = configLoad['Ab']  
    
    resultsTorpedo = getCapacityTorpedo(D1,D2,L1,L2,Hp,Su0,dSu,gamma,alphao,Np,rhows)
    
    print('******************  Torpedo Pile Result  *********************')

    print('Horizontal max. capacity,            ' , resultsTorpedo['Horizontal max.'], '[kN]')
    print('Vertical max. capacity,              ' , resultsTorpedo['Vertical max.'], '[kN]') 
    print('Unity check capacity,                ' , resultsTorpedo['UC'], '[-]') 

    print('**************************************************************')   
  


    
    
   
