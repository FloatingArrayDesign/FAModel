"""Drag embedment anchor capacity calculation functions in 'level 2' model,
currently set up for clay soils.
Lead author: Ahmed Radwan. 
"""

import numpy as np

def getCapacityDrag(Af, bm=0.073, En=1, Mu=0.1, Sum=0, N_eo=5.8, Theta_eo=58,
                    Xo_Lf=7.5, Zo=1, soil_type='clay', k=1.6):
    
    ''' 
    Calculate the inclined load capacity, measured in [kN], of the Drag Embedded Anchor in clay soil.
    The calculation is based on the soil properties, anchor geometry, and the mooring line properaties  
    
    Parameters
    ----------
    Af : float 
        The net area of the anchor's fluke [m2]
    bm : float 
        The mooring line diameter [m]   
    En : float 
        The multiplier number to be applied to chain bar diameter, if applicable (typical = 2.5, for wire line = 1)   
    Sum: float 
        The clay undrained shear strength at the mudline [kN/m2]
    k  : float 
        The gradient of the clay undrained shear strength with depth [kN/m2/m]   
    N_eo: float 
        The effective bearing capacity factor of the anchor     
    Theta_eo: float 
        The angle between the pad-eye load and the center of the fluke.    
    Xo_Lf: float 
        The intial ratio between the drag embedded distance to the fluke length.   
    Zo: float 
        The initial embedded depth for the anchor. [m]
   
        
    Returns
    -------
    max_To : flaot 
        The inlcind force capacity at the seabed [kN]
    Va : floar 
        The corresponding anchor steel volume [m3]
    '''
    
    if not soil_type == 'clay':
        raise Exception('Only clay soil type is supported for this anchor type so far.')
    
    # Anchor Geometry
    Wf_Lf = 1.25   # The aspect ratio of the plate width and the length
    Lf_tf = 7      # The aspect ratio of the plate length and the thickness
    Lf = np.sqrt(Af/Wf_Lf)   # The plate length, m
    Wf = Wf_Lf*Lf  # The plate width, m
    tf = Lf/Lf_tf  # The Plate Thickness, m

    # Shank Dimensions
    Ls_Lf = 1.25   # The aspect ratio of the shank length to the fluke length
    Ls = Ls_Lf*Lf  # The shank length, m
    ts_Lf = 0.14   # The aspect ratio of the shank thickness to the fluke length, m
    ts = ts_Lf*Lf  # The shank thickness, m.
    Ws = Ls/6
    ts_tf = ts/tf
    
    # Calculating the anchor steel volume
    Vf = Af*tf
    Vs = Ls*ts*Ws*2
    Va = Vf + Vs

    
    # The Anchor Initial Condition
    Suo = Sum + k*Zo  # The undrained shear strength at the initial embedded depth, kPa
    Nc = 10  # Bearing Capacity of the mooring line
    Tao = Af * N_eo * Suo
    Theta_a0 = np.sqrt((2*En*Nc*bm*Zo*(Sum + (k*Zo / 2)))/Tao)*180/np.pi

    Xo = Xo_Lf*Lf  # The initial drag distance, m

    Theta_f0 = Theta_eo - Theta_a0
    Rnto = 0.01
    dt = 0.1  # The tangential movement, m
    Su = Sum + k*Zo  # Initial value for Su


    # Specifying the initial condition calculations
    Theta_a = Theta_a0
    Theta_f = Theta_f0

    T_values = []  
    To_values = []  
    X_values = []  
    dZ_values = []  
    Theta_f_values = []  
    Zo_values = []  

    for _ in range(6000):
        dZ = (dt*np.sin(np.radians(Theta_f))) - (dt*Rnto*np.cos(np.radians(Theta_f)))
        Zo = round(Zo + dZ,3)
        Su = Sum + k*Zo
        dX = (dt*np.cos(np.radians(Theta_f))) + (dt*Rnto*np.sin(np.radians(Theta_f)))
        dX = dt * (1 + Rnto)*np.sin(np.radians(Theta_f))
        Xo = round(Xo + dX,3)
        F_dTheta_a = ((En*Nc*bm)/(N_eo*Af)) - ((k*(np.radians(Theta_a)**2))/(2*Su))
        dTheta_a = (F_dTheta_a*dZ/(Theta_a * np.pi / 180)) * (180 / np.pi)
        T =round(Af*N_eo*Su,2)
        Theta_a = Theta_a + dTheta_a
        To = T*np.exp(Mu*np.radians(Theta_a))  # holding capacity (at seabed surface)
        Theta_f = Theta_f - dTheta_a
        
        # Breaking the loop at when the fluke has been almost horizontal configuration
        Theta_ff= round(Theta_f, 0)
        if Theta_ff == 1:
           break       
       
        ts_tf = ts/tf
        dZ_values.append(dZ)
        Zo_values.append(Zo)
        T_values.append(T)
        To_values.append(To)
        X_values.append(Xo)
        Theta_f_values.append(Theta_f)  


    results = {}
    results['capacity'] = max(T_values)
    #results['vol'] = V
    results['embedment_depth'] = Zo
    results['drag_distance'] = Xo
    
    return results


if __name__ == '__main__':
 
    
    ''' Testing the function in one case of the net area of the anchor's fluke = 10 m2, 
    the aspect ratio of the plate width and its thickness, default is 40, all other parameters are default.'''
 
    results = getCapacityDrag(10, bm=0.073, En=1, Mu=0.1, Sum=0, N_eo=5.8, Theta_eo=58,
                        Xo_Lf=7.5, Zo=1, soil_type='clay', k=1.6)
    print('********************* One Case Test Result********************')

    print('Embedded length,              ' , results['embedment_depth'], '[m]') 
    print('Drag distance,                ' , results['drag_distance'], '[m]') 
    print('Load capacity,                ' , results['capacity'], '[kN]') 

    print('**************************************************************') 