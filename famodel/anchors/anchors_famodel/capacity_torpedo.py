
import yaml     
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from matplotlib import cm
from mpl_toolkits import mplot3d 
from capacity_load import getAnchorLoad
   
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
    # Xa = [5.82, 5.62, 5.04, 4.1, 2.91, 1.51, 0.00]
    # Ya = [0, 3.91846, 4.64614, 5.08714, 5.3706,	5.53709, 5.59723]
    # Xb = [11.6382, 11.2416,	10.0789, 8.22942, 5.81908, 3.01218, 7.12631E-16]
    # Yb = [0, 5.29435, 6.27754, 6.87338, 7.25637, 7.48133, 7.56259]
    
    fig, ax = plt.subplots()
    # ax.plot(Xa, Ya, color = 'r')
    # ax.plot(Xb, Yb, color = 'b')
    ax.plot(X, Y, color = 'g')
    ax.legend(['VHcapacity_soilA', 'VHcapacity_soilB', 'VHcapacity_soilC'])
    soilA = [(5.9,0.0),(5.9,1.6),(5.6,3.2),(4.7,4.7),(3.0,5.2),(1.4,5.1),(0.0,5.2)]
    soilB = [(11.7,0.0),(11.5,3.1),(10.6,6.1),(7.4,7.4),(4.5,7.7),(2.1,7.7),(0.0,7.8)]
    soilC = [(21.0,0.0),(20.9,5.6),(17.7,10.2),(11.11,11.11),(6.5,11.2),(3.1,11.4),(0.0,11.5)]
    xA = [x[0] for x in soilA]; yA = [x[1] for x in soilA]
    xB = [x[0] for x in soilB]; yB = [x[1] for x in soilB]
    xC = [x[0] for x in soilC]; yC = [x[1] for x in soilC]
    ax.scatter(xA, yA, s = 25, facecolors = 'none', edgecolors = 'r') # Soil A, ANSYS [Sousa, 2011]
    ax.scatter(xB, yB, s = 25, facecolors = 'none', edgecolors = 'b') # Soil B, ANSYS [Sousa, 2011]
    ax.scatter(xC, yC, s = 25, facecolors = 'none', edgecolors = 'g') # Soil C, ANSYS [Sousa, 2011]
       
    # Set labels and title
    plt.xlabel('Horizontal capacity, H [MN]')
    plt.ylabel('Vertical capacity, V [MN]')
    plt.suptitle('VH torpedo pile capacity envelope')
    #plt.axis([0,1.3*max(X[0],resultsLoad['H']),0,1.3*max(Y[-1],resultsLoad['V'])]) 
    plt.ylim(-2, 1.25*max(yC))
    plt.grid(True)
    plt.ioff()
    plt.show()

    fC = np.sqrt(X**2 + Y**2)     
    fAr = [np.sqrt((x[0])**2 + (x[1])**2) for x in soilA]
    fBr = [np.sqrt((x[0])**2 + (x[1])**2) for x in soilB]
    fCr = [np.sqrt((x[0])**2 + (x[1])**2) for x in soilC]  
       
    resultsTorpedo = {}
    resultsTorpedo['Horizontal max.'] = Hmax #Hmax[0]    # Capacity at specified loading angle
    resultsTorpedo['Vertical max.'] = Vmax               # Capacity at specified loading angle
   
    return resultsTorpedo