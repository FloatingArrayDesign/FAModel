# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:53:52 2024

@author: fmoreno
"""

import yaml      # Allow access to config file for user inputs
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def getAnchorLoad(Tm, thetam, zlug, Su0, dSu, nhu, Ab, Nc):
    
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
    Su0 : float 
        Undrained shear strength at the mudline (clay only), default is 2.39. [kPa]
    dSu : float 
        Undrained shear strength gradient (clay only), default is 1.41. [kPa/m]
    nhu : float 
        API- Factor [---], default is 0.5. (sand only)
    Ab : float
        Skin friction coefficient (clay only), default is 0.7.
    Nc : float
        Skin friction coefficient (clay only), default is 0.7.
    
    Returns
    -------
    Ta : flaot 
        Inclined load magnitude at the anchor lug [kN]
    thetaa : float 
        Inclined load angle at the anchor lug [deg]
    '''
    
    #xlug = 10;
    ez = (Su0*zlug**2/2 + dSu*zlug**3/3)/(Su0*zlug + dSu*zlug**2/2)
    ez_lug = ez/zlug
    Su_av_lug = Su0 + dSu*zlug/3
    zexQap = Ab*Nc*Su_av_lug
    
    def LoadTransfer(beta):
        return(np.sqrt((2*np.e**(nhu*(beta - thetam))*zexQap)/Tm + thetam**2) - beta)
    
    thetaa = fsolve(LoadTransfer,thetam); 
    thetaa = thetaa[0]; #print(thetaa)
    Ta = Tm/(np.e**(nhu*(thetaa - thetam)))
    
    H = Ta*np.cos(np.radians(thetaa)) 
    V = Ta*np.sin(np.radians(thetaa))
    
    # ax = plt.axes()
    # ax.arrow(xlug,0,Tm*np.cos(np.deg2rad(thetam))/100,Ta*np.sin(np.deg2rad(thetam))/100,head_width=5,head_length=10)
    # ax.arrow(0,-zlug,Ta*np.cos(np.deg2rad(thetaa))/100,Ta*np.sin(np.deg2rad(thetaa))/100,head_width=5,head_length=10)
    
    # Set labels and title
    # plt.xlabel('Horizontal capacity [kN]')
    # plt.ylabel('Vertical capacity [kN]')
    # plt.suptitle('VH suction pile capacity envelope')
    # plt.xlim(-1.2*xlug,1.25*Ta*np.cos(np.deg2rad(thetaa))/100); plt.ylim(-1.2*zlug,1.25*Ta*np.sin(np.deg2rad(thetaa))/100)
    # plt.grid(True)
    # plt.show()
    
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
     
