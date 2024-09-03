# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:53:52 2024

@author: fmoreno
"""

import yaml      # Allow access to config file for user inputs
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def getAnchorLoad(Tm, thetam, zlug, d, soil_type, Su0, k):
    
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
    d : float
        Chain diameter [m]
    
    Returns
    -------
    Ta : float 
        Inclined load magnitude at the anchor lug [kN]
    thetaa : float 
        Inclined load angle at the anchor lug [deg]
    '''
    
    # Setting bearing capacity values per soil type
    if soil_type == 'clay':
        Nc = 8.5; Ab=2.5; nhu=0.40
    elif soil_type == 'sand':
        Nc = 9; Ab=2.5; nhu=0.40
        
    thetam =  np.radians(thetam)
    print(thetam)
    
    if soil_type == 'clay':                           
        Su_av_lug = Su0*zlug + k*zlug**2/2
        zaQav = Ab*d*Nc*Su_av_lug
        
    elif soil_type == 'sand':
        gamma = 9
        zaQav = Ab*d*Nc*gamma*zlug**2/2
        
    def LoadTransfer(beta):
        return(2*zaQav*np.e**(nhu*(beta - thetam)) - Tm *(beta**2 - thetam**2))
    
    thetaa = fsolve(LoadTransfer, thetam) 
    thetaa = thetaa[0] 
    Ta = Tm/(np.e**(nhu*(thetaa - thetam)))
    
    H = FoSh*Ta*np.cos(thetaa) 
    V = FoSv*Ta*np.sin(thetaa)
       
    resultsLoad = {}
    resultsLoad['load'] = Ta                   # Load magnitude @ lug
    resultsLoad['angle'] = np.rad2deg(thetaa)  # Load angle @ lug
    resultsLoad['H'] = H                       # Horizontal component @ lug
    resultsLoad['V'] = V                       # Vertical component @ lug
    
    return resultsLoad

def getTransferLoad(Tm, thetam, zlug, line_type, d, soil_type, Su0, k):
    
    '''Calculate the transfer load from the mudline to the main padeye elevation using the DNV standards.
    The calculation is based on the mooring line properties, anchor geometry and the load from MoorPy and RAFT.  
    
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
    d : float
        Chain diameter [m]        
    
    Returns
    -------
    Ta : float 
        Inclined load magnitude at the anchor lug [kN]
    thetaa : float 
        Inclined load angle at the anchor lug [deg]
    '''
    
       
    # Setting bearing capacity values per soil type
    if soil_type == 'clay':
        Nc = 8.5
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
     
    resultsLoadDNV = {}
    resultsLoadDNV['load'] = Ta                   # Load magnitude @ lug
    resultsLoadDNV['angle'] = np.rad2deg(thetaa)  # Load angle @ lug
    resultsLoadDNV['H'] = H                       # Horizontal component @ lug
    resultsLoadDNV['V'] = V                       # Vertical component @ lug
    resultsLoadDNV['length'] = length_values      # Length of the embedded line
    
    # Plot of the line and extreme line tension
    drag_values = [-1*i for i in drag_values]
    depth_values = [-1*j for j in depth_values]
    fig, ax = plt.subplots(figsize=(20, 5)); n = 2000
    ax.scatter(drag_values[-1],depth_values[-1],color='g',zorder=5)
    ax.scatter(0,0,color='r',zorder=4)
    ax.arrow(0,0,Tm*np.cos(np.deg2rad(thetam))/n,Tm*np.sin(np.deg2rad(thetam))/n,
          head_width=0.25,head_length=0.5,color='r',zorder=3)
    ax.arrow(drag_values[-1],depth_values[-1],Ta*np.cos(thetaa)/n,Ta*np.sin(thetaa)/n,
          head_width=0.25,head_length=0.5,color='g',zorder=2)
    ax.plot(drag_values,depth_values,color='b',zorder=1)
     
    #Set labels and title
    plt.xlabel('Drag distance [m]')
    plt.ylabel('Embedded depth [m]')
    plt.suptitle('Inverse catenary profile in soil DNV')
    #plt.axis([0,1.3*max(X[0],H),0,1.3*max(Y[-1],V)]) 
    plt.grid(True)
    plt.show()
    
    return resultsLoadDNV
    
if __name__ == '__main__':
          
    ''' 
    Testing the function in one case
    '''
    '''
    # retrieves input data from a separate config file
    with open('LoadConfig.yml', 'r') as f:
        config = yaml.full_load(f)
    Tm = config['Tm']; thetam = config['thetam']
    zlug = config['zlug']; 
    Su0 = config['Su0']; dSu = config['dSu']
    nhu = config['nhu']; 
    Nc = config['Nc'];  Ab = config['Ab']  
    '''
    
    zlug = 14.31
    Tm = 21000
    thetam = 5
    
    #resultsLoad = getAnchorLoad(Tm, thetam, zlug)
    resultsLoadDNV = getAnchorLoadDNV(Tm, thetam, zlug)
    
    # print('******************* Load Transfer Resuls *********************')

    # print('Anchor load,                          ' , resultsLoad['load'], '[kN]') 
    # print('Anchor angle,                         ' , resultsLoad['angle'], '[deg]')
    # print('Horizontal load @ lug depth,          ' , resultsLoad['H'], '[kN]') 
    # print('Vertical load @ lug depth,            ' , resultsLoad['V'], '[kN]')

    # print('**************************************************************') 
     

          
    ''' 
    Testing the function in one case

    # retrieves input data from a separate config file
    with open('LoadConfig.yml', 'r') as f:
        config = yaml.full_load(f)
    Tm = config['Tm']; thetam = config['thetam']
    zlug = config['zlug']; 
    Su0 = config['Su0']; k = config['k']
    nhu = config['nhu']; 
    Nc = config['Nc']; Ab = config['Ab']  
    '''
    
    
    print('***************** Load Transfer Results DNV *******************')

    print('Anchor load,                          ' , resultsLoadDNV['load'], '[kN]') 
    print('Anchor angle,                         ' , resultsLoadDNV['angle'], '[deg]')
    print('Horizontal load @ lug depth,          ' , resultsLoadDNV['H'], '[kN]') 
    print('Vertical load @ lug depth,            ' , resultsLoadDNV['V'], '[kN]')

    print('**************************************************************') 

