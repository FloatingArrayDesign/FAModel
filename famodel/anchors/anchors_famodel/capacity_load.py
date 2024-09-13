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
    soil_type: string
        Select soil condition, 'clay' or 'sand'
    d : float
        Chain diameter [m]
    Su0 : float 
        Undrained shear strength at the mudline (clay) [kPa]
    k : float 
        Undrained shear strength gradient (clay) [kPa/m]
  
    Returns
    -------
    Ta : float 
        Inclined load magnitude at the anchor lug [kN]
    thetaa : float 
        Inclined load angle at the anchor lug [deg]
    '''
    
    # Setting bearing capacity values per soil type
    if soil_type == 'clay':
        Nc = 8.5; AB=2.5; nhu=0.40   # Nc    - Bearing capacity factor (9 and 14) DNV-RP-E301 
    elif soil_type == 'sand':        # Ab    - Effective unit bearing area (2.5 - 2.6 times chain dia)
        Nc = 9; AB=2.5; nhu=0.40     # nhu   - Friction coefficient between the mooring line and soil
        
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

def getTransferLoad(Tm, thetam, zlug, line_type, d, soil_type, Su0=None, k=None, w=None):

    
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
    soil_type = string
        Select soil condition, 'clay' or 'sand'
    d : float
        Chain diameter [m]    
    Su0 : float 
        Undrained shear strength at the mudline (clay only) [kPa]
    k : float 
        Undrained shear strength gradient (clay only) [kPa/m]
    w : float
        Mooring line unit weight [kg/m]
    
    Returns
    -------
    Ta : float 
        Inclined load magnitude at the anchor lug [kN]
    thetaa : float 
        Inclined load angle at the anchor lug [deg]
    '''
    
    deltas = 0.2
    
    # Setting bearing capacity values per soil type
    if soil_type == 'clay':
        Nc = 8.5
    elif soil_type == 'sand':
        Nc = 9
    
    # Include element weight in terms of d and match it with deltas
    if line_type == 'chain':   
        AS = 11.3*d; AB = 2.5*d; alpha = 0.5; W = w*deltas; 
    elif line_type == 'wire':
        AS = np.pi*d; AB = d; alpha = 0.3; W = w*deltas;
    
    T = Tm; theta = np.deg2rad(thetam);
    Su = Su0; drag = 0; depth = 0
         
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
    H = Ta*np.cos(thetaa); V = Ta*np.sin(thetaa) 
    length_values = deltas*len(drag_values)
     
    resultsLoad = {}
    resultsLoad['load'] = Ta                   # Load magnitude @ lug
    resultsLoad['angle'] = np.rad2deg(thetaa)  # Load angle @ lug
    resultsLoad['H'] = H                       # Horizontal component @ lug
    resultsLoad['V'] = V                       # Vertical component @ lug
    resultsLoad['length'] = length_values      # Length of the embedded line
    
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

    return resultsLoad

