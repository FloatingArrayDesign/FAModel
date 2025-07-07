# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:53:52 2024

@author: fmoreno
"""

import yaml      # Allow access to config file for user inputs
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def getAnchorLoad(Tm, thetam, zlug, d, soil_type, gamma, Su0, k):
    
    '''Calculate the inclined load capacity of a Suction Caisson Anchor in sand or clay.
    The calculation is based on the soil properties, anchor geometry, and the angle of inclined load.
    Offshore Geotechnical Engineering (Randolph , page 323)
    
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
        Nc = 8.5; Ab=2.5; nhu=0.40   # Nc    - Bearing capacity factor (9 and 14) DNV-RP-E301 
    elif soil_type == 'sand':        # Ab    - Effective unit bearing area (2.5 - 2.6 times chain dia)
        Nc = 9; Ab=2.5; nhu=0.35     # nhu   - Friction coefficient between the mooring line and soil
        
    thetam =  np.radians(thetam)
    
    if soil_type == 'clay':                           
        Su_av_lug = Su0*zlug + k*zlug**2/2
        zaQav = Ab*d*Nc*Su_av_lug
        
    elif soil_type == 'sand':
        zaQav = Ab*d*Nc*gamma*zlug**2/2
        
    def LoadTransfer(beta):
        return(2*zaQav*np.e**(nhu*(beta - thetam)) - Tm*(beta**2 - thetam**2))
    
    thetaa = fsolve(LoadTransfer, thetam) 
    thetaa = thetaa[0] 
    Ta = Tm/(np.e**(nhu*(thetaa - thetam)))
    
    H = Ta*np.cos(thetaa) 
    V = Ta*np.sin(thetaa)
       
    resultsLoad = {}
    resultsLoad['load'] = Ta                   # Load magnitude @ lug
    resultsLoad['angle'] = np.rad2deg(thetaa)  # Load angle @ lug
    resultsLoad['H'] = H                       # Horizontal component @ lug
    resultsLoad['V'] = V                       # Vertical component @ lug
    
    return resultsLoad

def getTransferLoad(Tm, thetam, zlug, line_type, d, soil_type, Su0=None, 
    k=None, gamma=None, phi= None, delta=None, w=None, plot=False):
    '''Calculate the transfer load from the mudline to the main padeye 
    elevation using the DNV standards. The calculation is based on the 
    mooring line properties, anchor geometry and the load from MoorPy and 
    RAFT.
    
    Parameters
    ----------
    Tm : float 
        Mooring line load at mudlevel [kN]
    thetam : float 
        Mooring line angle at mudlevel [deg]
    zlug : float 
        Embedded depth of the lug [m]
    line_type = string
        Select line type, 'chain' or 'wire'
    d : float
        Chain diameter [m]
    soil_type = string
        Select soil condition, 'clay' or 'sand'  
    Su0 : float 
        Undrained shear strength at the mudline (clay only) [Pa]
    k : float 
        Undrained shear strength gradient (clay only) [Pa/m]
    gamma: float
        Effective unit weight of the soil (sand only) [N/m3]
    phi : float
        Friction angle (sand only) [deg]
    delta: float
        Interface friction angle at soil-anchor line (sand only) [deg]
    w : float
        Mooring line unit weight [N/m]
    
    Returns
    -------
    Ta : float 
        Inclined load magnitude at the anchor lug [kN]
    thetaa : float 
        Inclined load angle at the anchor lug [deg]
    '''
    
    deltas = 0.2
      
    # Include element weight in terms of d and match it with deltas
    if line_type == 'chain':   
        Et = 10; En = 2.5;  W = w*deltas; 
    elif line_type == 'wire':
        Et = np.pi; En = 1; W = w*deltas;
    
    T = Tm; theta = np.deg2rad(thetam);
    Su = Su0; 
    drag = 0; depth = 0.1
         
    T_values = []; Su_values = [];
    drag_values = []; depth_values = []; 
    
    # Setting bearing capacity values per soil type
    if soil_type == 'clay':
        Nc = 8.5; alpha = 0.7;
    elif soil_type == 'sand':
        nhu = 0.5
        Nq = np.exp(np.pi*np.tan(np.deg2rad(phi)))*(np.tan(np.deg2rad(45 + phi/2)))**2
        # print(Nq)
        
    while (zlug - depth) >= 0:
        if soil_type =='clay':
            dtheta = (En*d*Nc*Su - W*np.cos(theta))/T*deltas
            dT = (Et*d*alpha*Su  + W*np.sin(theta))*deltas
            
        elif soil_type =='sand':
            dtheta  = (En*d*Nq*gamma*depth - W*np.cos(theta))/T*deltas
            dT = (Et*d*gamma*depth*np.tan(np.rad2deg(delta)) + W*np.sin(theta))*deltas 
        
        ddrag  = deltas*np.cos(theta)
        ddepth = deltas*np.sin(theta)
        theta += dtheta; T -= dT;
        
        drag += ddrag; depth += ddepth
        if Su:
            Su = Su0 + k*depth
                
        # Ensure consistency in load transfer
        if abs(Tm - T) > 0.75*Tm:  # More than 75% loss
            raise Exception(f"Warning: Load transfer is unrealistic. Initial load Tm = {Tm/1e6:.2f} MN and current load T = {T/1e6:.2f} MN differ by more than 75 %")
            break  # Exit the loop if load transfer is unrealistic
            
        # Check for excessive load angles
        if not (0 < np.rad2deg(theta) < 90):
            raise Exception(f"Warning: Load angle is unrealistic: {np.rad2deg(theta):.2f} deg")
            break  # Exit the loop if the angle becomes unreasonable   
            
        T_values.append(T); Su_values.append(Su)
        drag_values.append(drag); depth_values.append(depth); 
               
    Ta = T; thetaa = theta 
    # print(thetaa); print(Ta)
    H = Ta*np.cos(thetaa); V = Ta*np.sin(thetaa) 
    length_values = deltas*len(drag_values)
     
    resultsLoad = {}
    resultsLoad['diff'] = (Tm - Ta)/1e6        # Difference
    resultsLoad['load'] = Ta/1e6               # Load magnitude @ lug
    resultsLoad['angle'] = np.rad2deg(thetaa)  # Load angle @ lug
    resultsLoad['H'] = H                       # Horizontal component @ lug
    resultsLoad['V'] = V                       # Vertical component @ lug
    resultsLoad['length'] = length_values      # Length of the embedded line
    
    # Plot of the line and extreme line tension
    drag_values = [-1*i for i in drag_values]
    depth_values = [-1*j for j in depth_values]
    
    if plot:
        fig, ax = plt.subplots(figsize=(20, 5)); n = 2000000
        ax.scatter(drag_values[-1], depth_values[-1], color='g', zorder=5)
        ax.scatter(0, 0, color='r', zorder=4)
        ax.arrow(0, 0, Tm*np.cos(np.deg2rad(thetam))/n, Tm*np.sin(np.deg2rad(thetam))/n,
              head_width=0.25, head_length=0.5, color='r', zorder=3)
        ax.arrow(drag_values[-1], depth_values[-1], Ta*np.cos(thetaa)/n, Ta*np.sin(thetaa)/n,
              head_width=0.25, head_length=0.5, color='g',zorder=2)
        ax.plot(drag_values, depth_values,color='b', zorder=1)
         
        #Set labels and title
        plt.xlabel('Drag distance [m]')
        plt.ylabel('Embedded depth [m]')
        plt.suptitle('Inverse catenary profile in soil DNV')
        plt.grid(True)
    
    return resultsLoad


if __name__ == '__main__':
    
    Tm = 1.16e6 
    thetam = 0
    zlug = 10
    line_type ='chain'
    d = 0.160
    soil_type ='sand'
    Su0 = 2.4*1e3
    k = 1.41*1e3
    gamma = 9e3
    phi = 35
    delta = 27
    w = 4093

    resultsDNV = getTransferLoad(Tm, thetam, zlug, line_type, d, soil_type, Su0, k, gamma, phi, delta, w)
    #results = getAnchorLoad(Tm, thetam, zlug, d, soil_type, gamma, Su0, k)