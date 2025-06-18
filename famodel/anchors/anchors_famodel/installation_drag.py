
import yaml       # Allow access to config file for user inputs
import numpy as np
import matplotlib.pyplot as plt

def getInstallationDrag(Af, Lf, Ls, Lca, Lj, thetafs, bm, En, 
                    nhu, Su0, k, Ne, thetae0, z0,
                    Nn_max, Nt_max, Nm_max, m, n, p, q, plot=True): 
        
    ''' 
    Calculate the inclined load capacity, measured in [kN], of the Drag Embedded Anchor in clay soil.
    The calculation is based on the soil properties, anchor geometry, and the mooring line properaties  
    
    Parameters
    ----------
    Af : float 
        Net area of the anchor's fluke. [m2]
    Lf : float 
        Length of the anchor's fluke. [m]
    Ls : float 
        Length of the anchor's shank. [m]
    Lca : float 
        Perpendicular distance from the fluke to the shank. [m]
    Lj : float 
        Length of the anchor's fluke behind the shank. [m]
    bm : float 
        Mooring line diameter. [m]   
    En : float 
        Multiplier number to be applied to chain bar diameter, if applicable (typical = 2.5, for wire line = 1)   
    nhu : float
    Su0 : float 
        Clay undrained shear strength at the mudline. [kN/m2]
    k : float 
        Gradient of the clay undrained shear strength with depth. [kN/m2/m]   
    Ne : float 
        Effective bearing capacity factor of the anchor. [-]     
    theta_eo : float 
        Angle between the pad-eye load and the center of the fluke. [deg]     
    Z0 : float 
        Initial embedded depth for the anchor. [m]
    Nn_max, Nt_max, Nm_max : float 
        Bearing capacity factors. [-]
    n, m, p, q : float 
        Dimensionless coefficients. [-]
       
    Returns
    -------
    x : float 
        Drag embedded distance [m]
    z : float 
        Drag embedded depth [m]
    Qh : float 
        Holding capacity [kN]
    '''
        
    rhos = 7.850            # Steel density [t/m3]
    
    # Anchor Geometry
    Wf_Lf = 1.25            # Aspect ratio of the plate width and the length
    Lf_tf = 7               # Aspect ratio of the plate length and the thickness
    Lf = np.sqrt(Af/Wf_Lf)  # Plate length, m
    Wf = Wf_Lf*Lf           # Plate width, m
    tf = Lf/Lf_tf           # Plate thickness, m
    
    # Shank Dimensions
    Ls_Lf = 1.25            # Aspect ratio of the shank length to the fluke length
    Ls = Ls_Lf*Lf           # Shank length, m
    ts_Lf = 0.14            # Aspect ratio of the shank thickness to the fluke length, m
    ts = ts_Lf*Lf           # Shank thickness, m.
    Ws = Ls/6
    
    # Calculating the anchor steel weight
    Vf = Af*tf
    Vs = Ls*ts*Ws*2
    Va = round(Vf + Vs,1) 
    W = Va*rhos
    
    # The Anchor Initial Condition
    Su = Su0 + k*z0         # Undrained shear strength at the initial embedded depth, kPa
    Nc = 12                 # Bearing capacity of the mooring line
    Ta0 = Ne*Su*Af
    thetaa0 = np.sqrt((2*En*Nc*bm*z0*(Su0 + (k*z0/2)))/Ta0)
    thetae0 = np.deg2rad(thetae0)
    thetaf0 = thetae0 - thetaa0 
    
    x0 = 0; #X0_Lf*Lf       # Initial drag distance, m   
    Rnt0 = 0.01             # Uniform assumption of this value, -
    dt = 0.1                # Tangential movement, m
    
    # Specifying the initial condition calculations
    thetaa = thetaa0
    thetaf = thetaf0
    
    dz_values = []; z_values = []
    dx_values = []; x_values = []  
    Ta_values = []; To_values = [] 
    Su_values = [] 
    thetaa_values = []  
    thetaf_values = [] 
     
    z = z0; x = x0; Ta = Ta0
    xmax = 60; Tmax = 3000;
    
    for _ in range(3000):
        thetaaf = thetaf + thetaa
        thetaca = np.arctan(Ls*np.sin(np.deg2rad(thetafs))/(Lj + Ls*np.cos(np.deg2rad(thetafs)) - Lf/2))
        M = Ta*(np.cos(thetaaf)*np.sin(thetaca) - np.sin(thetaf)*np.cos(thetaca))*Lca
        Nn = Ta*np.sin(thetaaf); Nt = Ta*np.cos(thetaaf); Nm = M/(Su*Af*Lf)
        vn = (Nt_max/Nn_max)/(p*q/n)*(np.abs(Nn)/Nn_max)**(q-1)
        vt = ((np.abs(Nm)/Nm_max)**m + (np.abs(Nt)/Nt_max)**n)**(1/p-1)*(np.abs(Nt)/Nt_max)**(n-1) 
        Rnt = vn/vt
        # print(Rnt)
        dz = dt*np.sin(thetaf) - dt*Rnt*np.cos(thetaf)   
        dx = dt*np.cos(thetaf) + dt*Rnt*np.sin(thetaf)
        z += dz; x += dx
        dthetaa = (dz/thetaa)*(En*Nc*bm/Ne/Af - k*thetaa**2/2/Su)
        thetaa += dthetaa; thetaf -= dthetaa
        To = Ta*np.exp(nhu*thetaa) 
        Su = Su0 + k*z
        Ta = Ne*Su*Af
    
        dz_values.append(dz); z_values.append(z)   
        dx_values.append(dx); x_values.append(x)
        Ta_values.append(Ta); To_values.append(To)
        Su_values.append(Su)
        thetaa_values.append(thetaa) 
        thetaf_values.append(thetaf)  
        
        if thetaf < 0 or x > xmax or Ta > Tmax:
            break
    
    if plot:
        fig = plt.figure(figsize=(8, 6), dpi=100)
        plt.plot(x_values, z_values, color='red', linestyle='-')
        plt.gca().invert_yaxis()
        plt.xlabel('Drag distance (X), [m]')
        plt.ylabel('Embedded depth (z), [m]')
        plt.legend(['Drag Embedded Anchor (DEA) _ Clay Soil'], loc='upper left')
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.grid(True, which='major', linestyle='--')
        plt.show()
        
        fig = plt.figure(figsize=(8, 6), dpi=100)
        plt.plot(x_values, Ta_values, color='blue', linestyle='-')
        plt.xlabel('Drag distance (X), [m]')
        plt.ylabel('Maximum Tension Capacity (Ta), [kN]')
        plt.legend(['Drag Embedded Anchor (DEA) _ Clay Soil'], loc='upper left')
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.grid(True, which='major', linestyle='--')
        plt.show()
        
        fig = plt.figure(figsize=(8, 6), dpi=100)
        plt.plot(x_values, thetaf_values, color='green', linestyle='-')
        plt.xlabel('Drag distance (X), [m]')
        plt.ylabel('Fluke angle (deg), [kN]')
        plt.legend(['Drag Embedded Anchor (DEA) _ Clay Soil'], loc='upper left')
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.grid(True, which='major', linestyle='--')
        plt.show()
       
    resultsDrag = {}
    resultsDrag['capacity'] = max(Ta_values)
    resultsDrag['W'] = W
    resultsDrag['embedment_depth'] = z
    resultsDrag['drag_distance'] = x
    
    return resultsDrag
    
if __name__ == '__main__':
          
    ''' 
    Testing the function 
    '''
      
    # retrieves input data from a separate config file
    with open('DragConfig.yml', 'r') as f:
        configDrag = yaml.full_load(f)
    Af = configDrag['Af']; Lf = configDrag['Lf']; Ls = configDrag['Ls'];
    Lca = configDrag['Lf']; Lj = configDrag['Lj']; bm = configDrag['bm'];
    thetafs = configDrag['thetafs']; thetae0 = configDrag['thetae0'];
    x0_Lf = configDrag['x0_Lf']; En = configDrag['En']; Ne = configDrag['Ne']
    Su0 = configDrag['Su0']; k = configDrag['k']; nhu = configDrag['nhu'];
    Nn_max = configDrag['Nn_max']; Nt_max = configDrag['Nt_max']; Nm_max = configDrag['Nm_max'];
    m = configDrag['m']; n = configDrag['n']; p = configDrag['p']; q = configDrag['q'];
    z0 = configDrag['z0']; 
        
    resultsDrag = getInstallationDrag(Af, Lf, Ls, Lca, Lj, thetafs, bm,En, 
                                  nhu, Su0, k, Ne, thetae0, z0,
                                  Nn_max, Nt_max, Nm_max, m, n, p, q)  

    print('******************* Drag Anchor Result *****************')
    
    print('Embedded depth,               ' , resultsDrag['embedment_depth'], '[m]') 
    print('Drag distance,                ' , resultsDrag['drag_distance'], '[m]')
    print('Drag embedded anchor weight,  ' , resultsDrag['W'], '[t]')
    print('Load capacity,                ' , resultsDrag['capacity'], '[kN]') 
    
    print('********************************************************') 

