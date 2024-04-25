"""Suction caisson anchor capacity calculation functions in 'level 2' model, 
calculating inclined load capacity over the range of angles.
Lead author: Ahmed Radwan. 
"""

import numpy as np
import matplotlib.pyplot as plt
    

def getCapacitySuction(L,L_D_aspect=5 ,D_t_aspect=100, A_angle=60, F_angle= 2, 
                   soil_type='clay', gamma=0, Phi=30, So0= 2.39, k=1.41, Alpha=0.7, J=0.5):
    '''Calculate the inclined load capacity of a Suction Caisson Anchor in sand or clay.
    The calculation is based on the soil properties, anchor geometry, and the angle of inclined load.  
    
    Parameters
    ----------
    L : float 
        Suction anchor height [m]
    L_D_aspect : float 
        The aspect ratio of the anchor height and its outer diameter, default is 5.
    D_t_aspect : float 
        The aspect ratio of the anchor outer diameter and its thickness, default is 100.
    A_angle : float 
        The angle of the tension force at the padeye, default is 60. [deg] 
    F_angle : float 
        The adjustment factor of the load angle at the padeye, default is 2. [deg] 
    soil_type : string
        Specify 'sand' or 'clay'. This affects what other soil parameters are used.
    Phi : float
        The friction angle of the sand soil (sand only), default is 30. [deg] 
    gamma: float 
        The effective unit weight of the sand soil, default is 9 for sand, 4.7 for clay. [kN/m3] 
    J: float 
        The API- Factor [---], default is 0.5. (sand only)
    Su0 : float 
        The Undrained shear strength at the mudline (clay only), default is 2.39. [kPa]
    k : float 
        The Undrained shear strength gradient (clay only), default is 1.41. [kPa/m]
    Alpha : float
        The skin friction coefficient (clay only), default is 0.7. 
    
    Returns
    -------
    Ft_values : flaot 
        The inlcind force capacity at the padeye [kN]
    Anchor_vol : floar 
        The corresponding anchor steel volume [m3]
    '''
    
    # some constants
    specific_weight_steel = 78.5 # kN/m^3
    specific_weight_water = 9.81 # kN/m^3
    effective_weight_steel = specific_weight_steel - specific_weight_water
    
    # handle default gamma values
    if gamma==0:
        if soil_type=='clay':
            gamma = 4.7
        elif soil_type=='sand':
            gamma = 9
            
    # Suction Caisson Anchor Dimensions'
    D = round((L/L_D_aspect),2)                          # The diameter of the Anchor [m]
    t =  round(D/D_t_aspect, 3)                           # The thickness of the wall [m]
    Di = D - 2*t
    outer_r = D/2
    inner_r = (D - 2*t)/2
    height = L
    outer_vol = np.pi*outer_r**2*height
    inner_vol = np.pi*inner_r**2*height
    Anchor_vol = round(outer_vol - inner_vol,2)
    
    # ----- clay case -----
    if soil_type == 'clay':
    
        Li_percent = 0.67     # The pedeye location. 
        Li= L*Li_percent
        F_angle = 1           # NOTE: this is a hardcoded override of the input parameter <<<
    

        # Calculations related to the overburben pressure
        Lo_min = 0.1
        Lo_incr = 3
        N_inc = 34
        Lo = np.arange(Lo_min, Lo_incr*N_inc ,Lo_incr)
        F = abs(1 - (Li/Lo))
        n =100.00
        delta = L/n           # The increment of the soil profile. 
        inc_Z = L/100
        Z = np.arange(0.075, L, inc_Z)
        #fb=1.35
        
        # The calculations of the vertical capacity
        Su = So0 + k*Z
        Su_av = np.mean(Su)
        V = np.pi*D*L*Su_av*Alpha + 9*Su[-1]*np.pi*((D**2) /4)
        
        # The calculations of the horizontal capacity
        Np = 3 + gamma/Su + Z*(J/D)          # calculate Np element-wise
        Np = np.where(Np < 9, Np, 9)         # replace any values greater than 9 with 9
        Pult_D= Np * Su  
        VX_mat= np.zeros((len(Z), len(Lo)))
        for i in range(len(Lo)):
            VX_mat[:,i]=abs(1-(Z/Lo[i]))   
        dE_mat= np.zeros((len(Z), len(Lo)))  
        for i in range(len(Z)-1):
            for j in range(len(Lo)):
                dE_mat[i+1, j] = Pult_D[i+1]*delta*D*VX_mat[i+1, j] + dE_mat[i, j]        
        H = dE_mat[-1,:] / F
        Hm = min(H) #*fb

    
    # ----- sand case -----
    elif soil_type == 'sand':
       
        # The calculations of the Vertical capacity
        ktans = 0.5
        Zo = D/(4*ktans)
        Zi = Di/(4*ktans)
        Xo = L/Zo
        Xi = L/Zi
        Yo = np.exp(-Xo) - 1 + Xo
        Yi = np.exp(-Xi) - 1 + Xi
        V = gamma*(Zo**2)*Yo*ktans*np.pi*D + gamma*(Zi**2)*Yi*ktans*np.pi*Di 

        # The calculations of the horizontal capacity
        Nq = np.exp(np.pi*np.tan(np.radians(Phi)))*np.tan(np.radians(45 + (Phi/2)))**2
        Hm = 0.5*D*Nq*gamma*L**2

    else:
        raise Exception(f"Unsupported soil_type '{soil_type}'. Must be sand or clay.")
    
    
    # ----- add weight to total vertical capacity -----
    Vm = V + Anchor_vol*effective_weight_steel
        
        
    # ----- Interaction Curve (same for clay and sand) -----
    
    # Define the paramaters
    # Hm = Hmin   <<< removing because we can just use consistent variable names instead
    # Vm = Vm
    
    # Define the exponents
    a = 0.5 + L/D
    b = 4.5 + L/(3*D)
    
    # Define the function (for holding capacity?)
    def f(h, v):
        return (h/Hm)**a + (v/Vm)**b
    
    # Create a grid of horizontal and vertical load values
    h = np.linspace(0, Hm, 100)    # NOTE: what is the best choice for these values (clay was 100, sand was 2000) <<<
    v = np.linspace(0, Vm, 100)
    h_grid, v_grid = np.meshgrid(h, v)
    
    # Calculate the values of f on the grid
    f_grid = f(h_grid, v_grid)
    
    # Find the contour line where f < 1
    FB = min(0.99, (1**a) + (1**b))
    #contour_level = FB
    contour_lines = plt.contour(h_grid, v_grid, f_grid, levels=[FB], colors='b')
    #contour_lines = plt.contour(h_grid, v_grid, f_grid, colors='b')
    
    # Set labels and title
    plt.xlabel('Horizontal capacity [kN]')
    plt.ylabel('Vertical capacity [kN]')
    plt.suptitle('VH suction pile capacity envelope')
    plt.xlim(0,1.3*Hm); plt.ylim(0,1.3*Vm)
    
    # Extract the coordinates of the contour line
    contour_line = contour_lines.collections[0].get_paths()[0]
    h_values, v_values = contour_line.vertices[:, 0], contour_line.vertices[:, 1]
    
    # Estimating the inclined capacity
    angle = np.degrees(np.arctan(v_values/h_values))
    Ft  = np.sqrt(h_values**2 + v_values**2); print(Ft)
    #print(h_values); print(v_values)
    Ft1 = h_values/np.abs(np.cos(np.radians(angle)))
    Ft2 = v_values/np.abs(np.sin(np.radians(angle)))
    #print(Ft1); print(Ft2)
    
    # Extract the capacity values at the desired loading angle
    v_values_target  = np.interp(A_angle + F_angle, angle, v_values)
    h_values_target  = np.interp(A_angle + F_angle, angle, h_values)
    Ft_values_target = round(np.interp(A_angle + F_angle, angle, Ft),2)
    
    results = {}
    results['capacity'] = Ft_values_target  # capacity at specified loading angle
    results['capacity_force'] = Ft          # capacity envelope force magnitude
    results['capacity_angle'] = angle       # capacity envelope force angle
    results['vol'] = Anchor_vol
    results['L'] = L
    results['D'] = D
    results['t'] = t
    results['Anchor_vol'] = Anchor_vol
    
    return results
    
    

if __name__ == '__main__':
 
    
    ''' Testing the function in one case of the anchor with length = 15 m, 
    the aspect ratio of the anchor height and its outer diameter, default is 3, all other parameters are default.'''
 
    results = getCapacitySuction(15, L_D_aspect=5, soil_type='sand')
    print('********************* One Case Test Result********************')

    print('Anchor length,              ' , results['L'], '[m]') 
    print('Anchor diameter,            ' , results['D'], '[m]') 
    print('Anchor thickness,           ' , results['t'], '[m]')
    print('Anchor steel volume,        ' , results['Anchor_vol'], '[m3]') 
    print('Inclined load capacity,     ' , results['capacity'], '[kN]') 

    print('**************************************************************') 
