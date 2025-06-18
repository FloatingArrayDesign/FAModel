
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def py_Matlock(z, D, gamma, Su, sigma_v_eff, z0=None, return_curve=False):
   ''' Generate Matlock (1970) p–y curve at a given depth in clay.

    Parameters
    ----------
    z : float
        Depth relative to pile head (m)
    D : float
        Pile diameter (m)
    f_Su : function
        Undrained shear strength (Pa)
    f_sigma_v_eff : function
        Effective vertical stress (Pa)
    f_gamma : function
        Effective unit weight (kN/m³)
    z0 : float, optional
        Mudline depth (m). If provided, disables resistance above this level
    return_curve : bool
        Resulting p–y curve if True

    Returns
    -------
    f : interp1d
        Interpolation function for p–y relationship (N/m vs m)
    '''
    
    # Su = f_Su(z)
    # sigma_v_eff = f_sigma_v_eff(z)
    # gamma = f_gamma(z)
   
   # Strain at half the strength as defined by Matlock (1970). 
   # Typically ranges from 0.005 (stiff clay) to 0.02 (soft clay).
   epsilon_50 = 0.01          
   J = 0.5

   # No soil resistance above mudline
   if z0 is not None and z < z0:
       return lambda y_val: np.zeros_like(y_val)

   # Calculate ultimate resistance and shape parameters
   Nc = 3.0 + sigma_v_eff/Su + J*z/D
   Nc = min(Nc, 9.0)
   z_cr = 6.0 *D/(gamma*D/Su + J)

   p_ult = Su*Nc*D
   y_50 = 2.5*epsilon_50*D

   # Normalized lateral displacements
   N = 20
   # Y = np.concatenate((-np.logspace(3, -4, N), [0], np.logspace(-4, 3, N)))
   Y = np.linspace(-5, 5, 200)  # Normalized displacement in ±5 units
   P = 0.5 * np.sign(Y)*np.abs(Y)**(1.0/3.0)
   P = np.clip(P, -1.0, 1.0)

   # Un-normallized p-y curves
   y = Y*y_50
   p = P*p_ult

   f = interp1d(y, p, kind='linear', bounds_error=False, fill_value=0.0)   

   return (f, (y, p)) if return_curve else f 
       
def py_API(z, D, phi, sigma_v_eff, Dr, z0=None, return_curve=False):   
    ''' Generate API RP2A (1993) p–y curve at a given depth in sand.

    Parameters
    ----------
    z : float
        Depth relative to pile head (m)
    D : float
        Pile diameter (m)
    f_phi : function
        Friction angle (deg)
    f_sigma_v_eff : function
        Effective vertical stress (Pa)
    f_Dr : function
        Relative density (-)
    z0 : float, optional
        Mudline depth (m). If provided, disables resistance above this level
    plot : bool
        Resulting p–y curve if True

    Returns
    -------
    f : interp1d
        Interpolation function for p–y relationship (N/m vs m)
    '''
    
    # phi = f_phi(z)
    # sigma_v_eff = f_sigma_v_eff(z)
    # Dr = f_Dr(z)
         
    # Interpolate coefficients depending on the effective friction angle
    phi_ref = [  20,   25,   30,   35,   40]
    C1_ref  = [0.80, 1.25, 1.90, 3.00, 4.50]
    C2_ref  = [1.60, 2.10, 2.60, 3.40, 4.30]
    C3_ref  = [  10,   15,   30,   55,  105]
    
    C1 = np.interp(phi, phi_ref, C1_ref)
    C2 = np.interp(phi, phi_ref, C2_ref)
    C3 = np.interp(phi, phi_ref, C3_ref)
    
    # Disable p–y curve above mudline
    if z0 is not None and z < z0:
        return lambda y_val: np.zeros_like(y_val)
    
    # Compute ultimate lateral resistance
    p_ult = min(C1*z + C2*D, C3*D) * sigma_v_eff
    
    # Compute initial stiffness k (kN/m3 → N/m3)
    k = (54.6*Dr**2 + 0.8*Dr + 1.8)*1e3
    
    # Normalized displacement range
    N = 20
    y = np.concatenate((-np.logspace(3, -4, N), [0], np.logspace(-4, 3, N)))
    y = np.linspace(-0.05*D, 0.05*D, 200)  # ±5 cm for D = 1 m
    
    # Shape coefficient A
    A = max(3 - 0.8*z/D, 0.9)
    
    # Apply API p–y formulation
    ε = 1e-6  # prevent division by zero
    p = A*p_ult*np.tanh(k*z*y/(A*p_ult + ε))
    
    f = interp1d(y, p, kind='linear', bounds_error=False, fill_value=0.0)
    
    return (f, (y, p)) if return_curve else f  

def py_Reese(z, D, UCS, Em, z0=None, return_curve=False):
    ''' Generate Reese (1997) p–y curve at a given depth in weak rock.

    Parameters
    ----------
    z : float
        Depth relative to pile head (m)
    D : float
        Pile diameter (m)
    f_UCS : function
        Unconfined compressive strength UCS(z) (Pa)
    f_Em : function
        Young's modulus Em(z) (Pa)
    z0 : float, optional
        Mudline depth (m). If provided, disables resistance above this level
    plot : bool
        Resulting p–y curve if True

    Returns
    -------
    f : interp1d
        Interpolation function for p–y relationship (N/m vs m)
    '''
    
    # UCS = f_UCS(z)
    # Em = f_Em(z)
  
    RQD = 52                     # Assumed fair rock quality (moderately weathered rocks) 
    Dref = 0.305;                # Reference diamter (m)
    nhu = 0.3; E = 200e9
    t = (6.35 + D*20)/1e3        # Pile wall thickness (m), API RP2A-WSD
    I  = np.pi*(D**4 - (D - 2*t)**4)/64.0
    EI = E*I
    alpha = -0.00667*RQD + 1
    krm = 0.0005

    if z < z0:
        # Above mudline, no resistance
        p_ur = 0
    else:
        if z < 3*D:
            p_ur = alpha*UCS*D*(1 + 1.4*z/D)
            #kir = (100 +400*z/(3*D))
        else:
            p_ur = 5.2*alpha*UCS*D
            #kir = 500
        
    kir = (D/Dref)*2**(-2*nhu)*(EI/(Em*D**4))**0.284
    Kir = kir*Em     
    y_rm = krm*D
    y_a = (p_ur/(2*y_rm**0.25*Kir))**1.333    
   
    # Normalized lateral displacement
    N = 20
    y = np.concatenate((-np.logspace(4,-3,N),[0],np.logspace(-3,4,N)))
    y = np.linspace(-0.02*D, 0.02*D, 200)  # ±2 cm
    
    p = []
    for val in y:
        if abs(val) < y_a:
            p_val = np.sign(val)*Kir*val
        else:
            p_val = np.sign(val)*min((p_ur/2)*(abs(val)/y_rm)**0.25, p_ur)
        p.append(p_val)                                  
    
    f = interp1d(y, p, kind='linear', bounds_error=False, fill_value=0.0)
    
    return (f, (y, p)) if return_curve else f      

def py_Lovera(z, D, UCS, Em, zlug, z0, delta_grout=0.075, E_grout=20e9, delta_crushed=0.025, return_curve=False):
    ''' Generate Lovera (2019) p–y curve at a given depth for layered rock interfaces.

    Parameters
    ----------
    z : float
        Depth relative to pile head (m)
    D : float
        Pile diameter (m)
    f_UCS : function
        Unconfined compressive strength UCS(z) (Pa)
    f_Em : function
        Young's modulus Em(z) (Pa)
    zlug : float
        Load eccentricity (m)
    z0 : float
        Mudline depth (m). If provided, disables resistance above this level
    delta_grout : float, optional
        Grout annulus thickness (m) (default value: delta_grout=0.075)
    E_grout : float, optional
        Grout elastic modulus (Pa) (default value: E_grout=20e9)
    delta_crushed : float, optional
        Crushed rock annulus thickness (m) (default value: delta_crushed=0.025)
    plot : bool
        Resulting p–y curve if True

    Returns
    -------
    f : interp1d
        Interpolation function for p–y relationship (N/m vs m)
    '''
 
    if z < z0:
        return lambda y: np.zeros_like(y)

    # Retrieve elastic modulus at depth
    # Em = f_Em(z)
    nu = 0.3  # Typical Poisson's ratio for rock
    G_rock = Em/(2*(1 + nu))
    k_rock = 4*G_rock
    
    # Set E_crushed as 25% of intact rock modulus if not given
    E_crushed = 0.25*Em

    # Compute total stiffness from linear components
    k_eq = 1.0/(0.4*delta_grout/E_grout + delta_crushed/E_crushed + 1.0/k_rock)
    # print(f"z = {z:.2f}, Em = {Em:.2e}, k_eq = {k_eq:.2e} N/m")

    y = np.linspace(-0.03*D, 0.03*D, 200)
    p = k_eq*y
    
    f = interp1d(y, p, fill_value="extrapolate")

    return (f, (y, p)) if return_curve else f 

