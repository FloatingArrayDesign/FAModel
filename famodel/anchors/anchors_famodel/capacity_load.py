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

def getTransferLoad(Tm, thetam, zlug, line_type, d, soil_type, Su0=None, k=None, qc0=None, kc=None, w=None):

    
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
        Mooring line unit weight [kN/m]
    
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
        AS = 11.3*d; AB = 2.5*d; alpha = 0.7; W = w*deltas; 
    elif line_type == 'wire':
        AS = np.pi*d; AB = d; alpha = 0.3; W = w*deltas;
    
    T = Tm; theta = np.deg2rad(thetam);
    Su = Su0; qc = qc0;
    drag = 0; depth = 0
         
    T_values = []; Su_values = [];
    drag_values = []; depth_values = []; 
    
    # Setting bearing capacity values per soil type
    if soil_type == 'clay':
        Nc = 8.5
    elif soil_type == 'sand':
        EWB = 2.2; beta = 0.625; nhu = 0.46
        
    while (zlug - depth) >= 0:
        if soil_type =='clay':
            dtheta = (Nc*Su*AB - W*np.cos(theta))/T*deltas
            dT = -(alpha*Su*AS - W*np.sin(theta))*deltas
        if soil_type =='sand':
            dtheta = (EWB*d*beta*qc*AB - W*np.cos(theta))/T*deltas
            dT = -(nhu*EWB*d*beta*qc*AS - W*np.sin(theta))*deltas 
        ddrag = deltas*np.cos(theta)
        ddepth = deltas*np.sin(theta)
        theta += dtheta; T += dT; 
        drag += ddrag; depth += ddepth
        Su = Su0 + k*depth
        qc = qc0 + kc*depth
           
        T_values.append(T); Su_values.append(Su)
        drag_values.append(drag); depth_values.append(depth); 
               
    Ta = T; thetaa = theta 
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
    fig, ax = plt.subplots(figsize=(20, 5)); n = 2000000
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
    plt.grid(True)
    plt.show()   

    return resultsLoad

if __name__ == '__main__':
    
    Tm = 1.16e6  
    thetam = 0
    zlug = 10
    line_type ='chain'
    d = 0.160
    soil_type ='clay'
    Su0 = 2.4*1e3
    k = 1.41*1e3
    qc0 = 1.3e6
    kc = 1.5*1e6
    w = 4093

    resultsDNV = getTransferLoad(Tm, thetam, zlug, line_type, d, soil_type, Su0, k, qc0, kc, w)
    #results = getAnchorLoad(Tm, thetam, zlug, d, soil_type, gamma, Su0, k)