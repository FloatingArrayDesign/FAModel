# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:53:52 2024

@author: fmoreno
"""

import yaml      # Allow access to config file for user inputs
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def getAnchorLoad(Tm, thetam, zlug, safety_factor='yes', Su0=10.0, k=2.0, nhu=0.5, Nc=14, Ab=2.6):
    
    '''Calculate the inclined load capacity of a Suction Caisson Anchor in sand or clay.
    The calculation is based on the soil properties, anchor geometry, and the angle of inclined load.  
    
    Parameters
    ----------
    Tm : float 
        Mooring line load at mudlevel [kN]
    thetam : float 
        Mooring line angle at mudlevel [deg]
    zlug : float 
        Embedded depth of the lug [m]
    safety_factor: float
        Specify 'yes' or 'no'. This affects load and resistance factors.
    Su0 : float 
        Undrained shear strength at the mudline (clay) [kPa]
    k : float 
        Undrained shear strength gradient (clay) [kPa/m]
    nhu : float 
        Friction coefficient between the mooring line and soil
    Ab : float
        Effective unit bearing area (2.5 - 2.6 times Chain Dia)
    Nc : float
        Bearing capacity factor (9 and 14) DNV-RP-E301 
    
    Returns
    -------
    Ta : float 
        Inclined load magnitude at the anchor lug [kN]
    thetaa : float 
        Inclined load angle at the anchor lug [deg]
    '''
    
    # Setting default safety factors
    # API RP 2SK permanent intact condition lateral/axial safety factors. [-]
    if safety_factor == 'yes':
        FoSh = 1.6
        FoSv = 2.0
    elif safety_factor == 'no':
        FoSh = 1.0
        FoSv = 1.0
        
    ez = (Su0*zlug**2/2 + k*zlug**3/3)/(Su0*zlug + k*zlug**2/2)
    Su_av_lug = Su0 + k*zlug/3
    zexQap = Ab*Nc*Su_av_lug
    
    def LoadTransfer(beta):
        return(np.sqrt((2*np.e**(nhu*(beta - thetam))*zexQap)/Tm + thetam**2) - beta)
    
    thetaa = fsolve(LoadTransfer,thetam); 
    thetaa = thetaa[0]; #print(thetaa)
    Ta = Tm/(np.e**(nhu*(thetaa - thetam)))
    
    H = FoSh*Ta*np.cos(np.radians(thetaa)) 
    V = FoSv*Ta*np.sin(np.radians(thetaa))
       
    resultsLoad = {}
    resultsLoad['load'] = Ta                   # Capacity at specified loading angle
    resultsLoad['angle'] = thetaa              # Horizontal load @ lug
    resultsLoad['H'] = H                       # Capacity at specified loading angle
    resultsLoad['V'] = V                       # Horizontal load @ lug
    
    return resultsLoad
    
if __name__ == '__main__':
          
    ''' 
    Testing the function in one case
    '''
  
    # retrieves input data from a separate config file
    with open('LoadConfig.yml', 'r') as f:
        config = yaml.full_load(f)
    Tm = config['Tm']; thetam = config['thetam']
    zlug = config['zlug']; 
    Su0 = config['Su0']; k = config['k']
    nhu = config['nhu']; 
    Nc = config['Nc']; Ab = config['Ab']  
    
    resultsLoad = getAnchorLoad(Tm=5000, thetam=25, zlug=14)
    
    print('******************  Suction Pile Result  *********************')

    print('Anchor factored load,                 ' , resultsLoad['load'], '[kN]') 
    print('Anchor factored angle,                ' , resultsLoad['angle'], '[deg]')
    print('Horizontal factored load @ lug depth, ' , resultsLoad['H'], '[kN]') 
    print('Horizontal factored load @ lug depth, ' , resultsLoad['V'], '[kN]')

    print('**************************************************************') 
     

def getAnchorLoadDNV(Tm, thetam, zlug, safety_factor='yes', line_type='chain', soil_type='clay', Su0=10.0, k=2.0, d=0.15):
    
    '''Calculate the inclined load capacity of a Suction Caisson Anchor in sand or clay.
    The calculation is based on the soil properties, anchor geometry, and the angle of inclined load.  
    
    Parameters
    ----------
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
    Nc : float
        Bearing capacity factor (clay and sand)
    
    Returns
    -------
    Ta : float 
        Inclined load magnitude at the anchor lug [kN]
    thetaa : float 
        Inclined load angle at the anchor lug [deg]
    '''
    
    # Setting default safety factors
    # API RP 2SK permanent intact condition lateral/axial safety factors. [-]
    if safety_factor == 'yes':
        FoSh = 1.6
        FoSv = 2.0
    elif safety_factor == 'no':
        FoSh = 1.0
        FoSv = 1.0
        
    # Setting default gamma values per soil type
    # Effective unit weight of the sand soil, default is 9 for sand, 4.7 for clay. [kN/m3]
    if soil_type == 'clay':
        Nc = 11.5
    elif soil_type == 'sand':
        Nc = 9
    
    # Include element weight in terms of d and match it with deltas
    if line_type == 'chain':   
        AS = 11.3*d; AB = 2.5*d; alpha = 0.5; W = 15; 
    elif line_type == 'wire':
        AS = np.pi*d; AB = d; alpha = 0.3; W = 80;
    
    T = Tm; theta = np.deg2rad(thetam);
    Su = Su0; drag = 0; depth = 0
     
    deltas = 0.2
    
    T_values = []; Su_values = [];
    drag_values = []; depth_values = []; 
    
    while (zlug - depth) >= 0:
        dtheta = (Nc*Su*AB - W*np.cos(theta))/T*deltas
        dT = -(alpha*Su*AS - W*np.sin(theta))*deltas
        ddrag = deltas*np.cos(theta)
        ddepth = deltas*np.sin(theta)
        theta += dtheta; T += dT; 
        drag += ddrag; depth += ddepth
        Su = Su0 + k*depth
           
        T_values.append(T); Su_values.append(Su)
        drag_values.append(drag); depth_values.append(depth); 
               
    Ta = T; thetaa = theta 
    H = FoSh*Ta*np.cos(thetaa); V = FoSv*Ta*np.sin(thetaa) 
    length_values = deltas*len(drag_values)
     
    resultsLoad = {}
    resultsLoad['load'] = Ta                   # Capacity at specified loading angle
    resultsLoad['angle'] = np.rad2deg(thetaa)  # Horizontal load @ lug
    resultsLoad['H'] = H                       # Capacity at specified loading angle
    resultsLoad['V'] = V                       # Horizontal load @ lug
    resultsLoad['length'] = length_values      # Length of the embedded line
    
    # Plot of the line and extreme line tension
    drag_values = [-1*i for i in drag_values]
    depth_values = [-1*j for j in depth_values]
    ax = plt.axes()
    ax.scatter(drag_values[-1],depth_values[-1],color='g',zorder=5)
    ax.scatter(0,0,color='r',zorder=4)
    ax.arrow(0,0,Tm*np.cos(np.deg2rad(thetam))/500,Tm*np.sin(np.deg2rad(thetam))/500,
          head_width=0.25,head_length=0.5,color='r',zorder=3)
    ax.arrow(drag_values[-1],depth_values[-1],Ta*np.cos(thetaa)/500,Ta*np.sin(thetaa)/500,
          head_width=0.25,head_length=0.5,color='g',zorder=2)
    ax.plot(drag_values,depth_values,color='b',zorder=1)
     
    #Set labels and title
    plt.xlabel('Drag distance [m]')
    plt.ylabel('Embedded depth [m]')
    plt.suptitle('Inverse catenary profile in soil')
    plt.grid(True)
    plt.show()
    
    return resultsLoad
    
if __name__ == '__main__':
          
    ''' 
    Testing the function in one case
    '''
  
    # retrieves input data from a separate config file
    with open('LoadConfig.yml', 'r') as f:
        config = yaml.full_load(f)
    Tm = config['Tm']; thetam = config['thetam']
    zlug = config['zlug']; 
    Su0 = config['Su0']; dSu = config['dSu']
    nhu = config['nhu']; 
    Nc = config['Nc'];  Ab = config['Ab']  
    
    resultsLoad = getAnchorLoad(Tm, thetam, zlug, Su0, dSu, nhu, Ab, Nc)
    
    print('******************  Suction Pile Result  *********************')

    print('Anchor factored load,                 ' , resultsLoad['load'], '[kN]') 
    print('Anchor factored angle,                ' , resultsLoad['angle'], '[deg]')
    print('Horizontal factored load @ lug depth, ' , resultsLoad['H'], '[kN]') 
    print('Horizontal factored load @ lug depth, ' , resultsLoad['V'], '[kN]')

    print('**************************************************************') 
     
