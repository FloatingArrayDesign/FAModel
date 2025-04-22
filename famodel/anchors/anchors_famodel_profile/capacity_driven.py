
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import linalg
from capacity_soils import clay_profile, sand_profile, rock_profile
from capacity_pycurves import py_Matlock, py_API, py_Reese
from capacity_plots import plot_pile

def getCapacityDriven(profile, soil_type, L, D, zlug, H, V, plot=True):   
    '''Models a laterally loaded pile using the p-y method. The solution for
    lateral displacements is obtained by solving the 4th order ODE, EI*d4y/dz4
    EI*d4y/dz4 - V*d2y/dz2 + ky = 0 using the finite difference method.
    EI*d4y/dz4 - V*d2y/dz2 + K*z*dy/dz + ky = 0 using the finite difference method.
    
    Assumes that EI remains constant with respect to curvature i.e. pile
    material remains in the elastic region.

    Parameters
    ----------
    profile : array
        Soil profile as a 2D array: (z, parameters)
            Clay: (z (m), Su (kPa), gamma (kN/m³))
            Sand: (z (m), phi (deg), gamma (kN/m³), Dr (%))
            Rock: (z (m), UCS (MPa), Em (MPa))
    soil_type : string
        Select soil condition: 'clay', 'sand', or '(weak) rock'
    L : float
        Pile length (m)
    D : float
        Pile diameter (m)
    zlug : float
        Depth of padeye from pile head (m)
    H : float
        Horizontal load applied at padeye (N)
    V : float
        Vertical load applied at padeye (N)
    plot : bool
        Plot the p-y curve and the deflection pile condition if True

    Returns
    -------
    y : array
        Lateral displacement at each node (real nodes only)
    z : array
        Node depth positions corresponding to y (m)
    resultsDriven : dict
        Dictionary containing displacements, moment capacity, hinge state and vertical capacity
    '''

    n = 50; iterations = 10; loc = 2
    nhuc = 1; nhu = 0.3         # Resistance factor (-)
    delta_r = 0.08              # Mean roughness height (m)
    
    t = (6.35 + D*20)/1e3       # Pile wall thickness (m), API RP2A-WSD
    E = 200e9                   # Elastic modulus of pile material (Pa)
    fy = 350e6                  # Yield strength of pile material (Pa)
    rhows = 66.90e3             # Submerged steel specific weight (N/m3)
    rhow = 10e3                 # Water specific weight (N/m3) 
    
    # Pile geometry
    I = (np.pi/64.0)*(D**4 - (D - 2*t)**4)
    EI = E*I
    h  = L/n                    # Element size
    N  = (n + 1) + 4            # (n+1) Real + 4 Imaginary nodes
    
    # Outer and inner surface of the pile skirt
    def PileSurface(Len, Dia):
        Sp = np.pi*Dia*Len
        return Sp    
    # Dry and wet mass of the pile    
    def PileWeight(Len, Dia, tw, rho):
        Wp = ((np.pi/4)*(Dia**2 - (Dia - 2*tw)**2)*Len)*rho 
        return Wp 
    # Mass of the soil plug      
    def SoilWeight(Len, Dia, tw, gamma_soil): 
        Wsoil =(np.pi/4)*(Dia - 2*tw)**2*Len*gamma_soil
        return Wsoil

    # Array for displacements at nodes, including imaginary nodes.
    y = np.ones(N)*(0.01*D)     # An initial value of 0.01D was arbitrarily chosen

    # Initialize and assemble array/list of p-y curves at each real node
    z = np.zeros(N)
    py_funs  = []; PileShaft =[]
    k_secant = np.zeros(N)
    DQ = []

    for i in [0, 1]:            # Top two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0

    for i in range(2, n+3):    # Real nodes
        z[i] = (i - 2)*h      
        if soil_type == 'clay':
            z0, f_Su, f_sigma_v_eff, f_gamma, f_alpha = clay_profile(profile)
            if z[i] < z0:
                # No p-y curve above mudline
                py_funs.append(lambda y_val: np.zeros_like(y_val))
                k_secant[i] = 0.0
                PileShaft.append(0.0)
            else:
                Su = f_Su(z[i])
                sigma_v_eff = f_sigma_v_eff(z[i])
                gamma = f_gamma(z[i])
                alpha = f_alpha(z[i])
                py_funs.append(py_Matlock(z[i], D, zlug, f_Su, f_sigma_v_eff, f_gamma, z0=z0, plot=plot))
                Vo = np.pi*D*alpha*Su*z[i]**2
                PileShaft.append(Vo)
                k_val = py_funs[i](y[i])
                k_secant[i] = k_val/y[i] if y[i] != 0 else 0.0
            
        elif soil_type == 'sand':
            z0, f_phi, f_sigma_v_eff, f_gamma, f_Dr, f_delta = sand_profile(profile)
            if z[i] < z0:
                # No p-y curve above mudline
                py_funs.append(lambda y_val: np.zeros_like(y_val))
                k_secant[i] = 0.0
                PileShaft.append(0.0)
            else:
                phi = f_phi(z[i])
                sigma_v_eff = f_sigma_v_eff(z[i])
                gamma = f_gamma(z[i])
                delta = f_delta(z[i])
                py_funs.append(py_API(z[i], D, zlug, f_phi, f_sigma_v_eff, f_Dr, z0=z0, plot=plot))
                fs = delta * sigma_v_eff
                Vo = np.pi*D*fs*z[i]
                PileShaft.append(Vo)
                k_val = py_funs[i](y[i])
                k_secant[i] = k_val/y[i] if y[i] != 0 else 0.0
            
        elif soil_type in ['rock', 'weak_rock']:
            z0, f_UCS, f_Em = rock_profile(profile)            
            if z[i] < z0:
                # No p-y curve above mudline
                py_funs.append(lambda y_val: np.zeros_like(y_val))
                k_secant[i] = 0.0
                DQ.append(0.0)
            else:
                UCS = f_UCS(z[i])
                Em = f_Em(z[i])
                py_funs.append(py_Reese(z[i], D, zlug, f_UCS, f_Em, z0=z0, plot=plot))
                SCR = nhuc*Em/(UCS*(1 + nhu))*delta_r/D
                alpha = 0.36*SCR - 0.0005
                fs = alpha*UCS
                Dq = np.pi*D*fs*z[i]
                DQ.append(Dq)
                k_val = py_funs[i](y[i])
                k_secant[i] = k_val/y[i] if y[i] != 0 else 0.0
       
    for i in [n+3, n+4]:       # Bottom two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0

    # Compute individual contributions to total vertical load
    Wp    = PileWeight(L, D, t, rhows)                                           # Pile self-weight (wet)
    Wsoil = SoilWeight(L, D, t, gamma) if soil_type in ['clay', 'sand'] else 0.0 # Soil plug only in soil profiles
    Wshaft = PileShaft[-1] if soil_type in ['clay', 'sand'] else 0.0             # Shaft resistance for soils
    Wtip   = DQ[-1] if soil_type in ['rock', 'weak_rock'] else 0.0               # Tip resistance for rock
    
    # Compute total vertical capacity
    Vmax = Wp + Wsoil + Wshaft + Wtip
    
    for j in range(iterations):
        y, Mi, Mp, hinge_formed, hinge_location = fd_solver(n, N, h, D, t, fy, EI, H, V, zlug, z0, k_secant)
        for i in range(2, n+3):
            if callable(py_funs[i]):
                k_secant[i] = py_funs[i](y[i])/y[i]
                
    if plot:
        y1 = np.linspace(-2.*D, 2.*D, 500)
        plt.plot(y1, py_funs[loc](y1))
        plt.xlabel('y (m)'), plt.ylabel('p (N/m)')
        plt.grid(True)
        
        fig, ax = plt.subplots(figsize=(3, 5))
        y0 = np.zeros_like(z[2:-2])
        # Plot original vertical pile
        ax.plot(y0, z[2:-2], 'k', label='Original pile axis')
        # Plot deflected shape
        ax.plot(y[2:-2], z[2:-2], 'r', label='Deflected shape')
        # Padeye marker
        ax.plot(0, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')
        # Mudline marker 
        ax.axhline(z0, color='blue', linestyle='--', label=f'Mudline (z0 = {z0:.2f} m)')
        
        ax.set_xlabel('Lateral displacement (m)')
        ax.set_ylabel('Depth (m)')
        ax.set_xlim([-0.1*D, 0.1*D])
        ax.set_ylim([L + 5, -2])  # Downward is positive z 
        ax.grid(ls='--')
        ax.legend()           

    resultsDriven = {
        'Soil type': soil_type,
        'Weight pile': PileWeight(L, D, t, rhows + rhow),
        'Vertical max.': Vmax
    }
    
    if zlug <= 0:
    
        resultsDriven['Lateral displacement'] = y[2]
        resultsDriven['Rotational displacement'] = np.rad2deg((y[2] - y[3])/h)
    
    else:
        resultsDriven.update({
            'Lateral displacement': max(y),
            'Bending moment': Mi,
            'Plastic moment': Mp,
            'Plastic hinge': hinge_formed,
            'Hinge location': hinge_location
        })
    
    # Add model metadata
    if soil_type == 'clay':
        resultsDriven['p-y model'] = 'Matlock (1970)'
        resultsDriven['Soil profile'] = 'Su vs depth'
    elif soil_type == 'sand':
        resultsDriven['p-y model'] = 'API RP2A (1993)'
        resultsDriven['Soil profile'] = 'phi vs depth'
    elif soil_type in ['rock', 'weak_rock']:
        resultsDriven['p-y model'] = 'Reese (1997)'
        resultsDriven['Soil profile'] = 'UCS vs depth'

    
    return y[2:-2], z[2:-2], resultsDriven

def fd_solver(n, N, h, D, t, fy, EI, H, V, zlug, z0, k_secant):
    '''Solves the finite difference equations from 'py_analysis_1'. This function should be run iteratively for
    non-linear p-y curves by updating 'k_secant' using 'y'. A single iteration is sufficient if the p-y curves
    are linear.


    Parameters
    ----------
    n : int
        Number of elements
    N : int
        Total number of nodes (real + imaginary)
    h : float
        Element size (m)
    D : float
        Pile diameter (m)
    t : float
        Pile wall thickness (m)
    fy : float
        Yield strength of pile material (Pa)
    EI : float
        Flexural rigidity of the pile (Nm²)
    H : float
        Horizontal load at padeye (N)
    V : float
        Vertical load at padeye (N)
    zlug : float
        Depth of padeye from pile head (m)
    z0 : float
        Mudline depth from pile head (m)
    k_secant : array
        Secant stiffness from p-y curves at each node

    Returns
    -------
    y : array
        Lateral displacement at each node
    Mi : float
        Maximum internal bending moment (Nm)
    Mp : float
        Plastic moment capacity of the pile (Nm)
    hinge_formed : bool
        Whether plastic hinge is formed
    hinge_location : int
        Index of the node with hinge formation
    '''
     
    # Initialize and assemble matrix
    X = np.zeros((N, N))

    # (n+1) finite difference equations for (n+1) real nodes
    for i in range(0, n+1):
        X[i, i]   =  1.0
        X[i, i+1] = -4.0 + V*h**2/EI
        X[i, i+2] =  6.0 - 2*V*h**2/EI + k_secant[i+2]*h**4/EI
        X[i, i+3] = -4.0 + V*h**2/EI
        X[i, i+4] =  1.0

    # Curvature at pile head
    X[n+1, 1] =  1.0
    X[n+1, 2] = -2.0
    X[n+1, 3] =  1.0

    # Shear at pile head
    X[n+2, 0] = -1.0
    X[n+2, 1] =  2.0 - V*h**2/EI
    X[n+2, 2] =  0.0
    X[n+2, 3] = -2.0 + V*h**2/EI
    X[n+2, 4] =  1.0

    # Curvature at pile tip
    X[n+3, -2] =  1.0
    X[n+3, -3] = -2.0
    X[n+3, -4] =  1.0

    # Shear at pile tip
    X[n+4, -1] =  1.0
    X[n+4, -2] = -2.0 + V*h**2/EI
    X[n+4, -3] =  0.0
    X[n+4, -4] =  2.0 - V*h**2/EI
    X[n+4, -5] = -1.0

    # Initialize vector q
    q = np.zeros(N)
    
    # Always apply shear
    # Index of the node where the horizontal load is applied (padeye)
    zlug_index = int(zlug/h)
    q[zlug_index] = 2*H*h**3
                 
    y = linalg.solve(EI*X, q)

    # Compute the plastic moment capacity Mp
    Zp = (1/6)*(D**3 - (D - 2*t)**3)  # Plastic section modulus for hollow pile (m3)
    Mp = Zp*fy                        # Plastic moment capacity (N/m)

    # Check for plastic hinge formation
    Mi, Mp, hinge_formed, hinge_location = plastic_hinge(y, h, EI, Mp)
    
    return y, Mi, Mp, hinge_formed, hinge_location
   
def plastic_hinge(y, h, EI, Mp):
    '''Check for plastic hinge formation along the pile.
    
    Parameters
    ----------
    y : array
        Lateral displacements at each node
    h : float
        Element size (m)
    EI : float
        Flexural rigidity of the pile (Nm²)
    Mp : float
        Plastic moment capacity (Nm)

    Returns
    -------
    Mi_max : float
        Maximum internal moment along the pile (Nm)
    Mp : float
        Plastic moment capacity (Nm)
    hinge_formed : bool
        True if plastic hinge is formed
    hinge_location : int
        Node index where hinge forms (if any)
    '''
    
    hinge_formed = False
    hinge_location = -1
    Mi_all = []
        
    # Loop through each internal node and compute the bending moment
    for i in range(1, len(y) - 1):
        Mi = EI * (y[i+1] - 2*y[i] + y[i-1])/h**2
        Mi_all.append(Mi)

        if abs(Mi) >= Mp and not hinge_formed:
            hinge_formed = True
            hinge_location = i

    Mi_max = max(Mi_all, key=abs) if Mi_all else 0.0

    return Mi_max, Mp, hinge_formed, hinge_location

if __name__ == '__main__':
    
    profile_clay = np.array([
        [ 1.0, 600, 8],     
        [ 6.0, 200, 8],
        [15.0, 400, 8],
        [30.0, 600, 9]
    ])
    
    profile_sand = np.array([
        [ 2.0, 28,  8, 75],
        [10.0, 34,  9, 75],
        [15.0, 36, 10, 75],
        [40.0, 45,  9, 85]
    ])
    
    profile_rock = np.array([
        [ 2.0, 0.5, 1e3],
        [ 5.0, 2.0, 2e4],
        [30.0, 1.0, 5e4]
    ])

    D = 2.0           # Diameter (m)
    L = 15.0          # Length (m)
    zlug = 1          # Padeye depth (m)
    
    # === CLAY ===
    # y_clay, z_clay, results_clay = getCapacityDriven(profile_clay, 'clay', L, D, zlug, H=5.0e6,  V=1.5e5, plot=True)

    # plot_pile(profile_clay, 'clay', y_clay, z_clay, D, L, profile_clay[0, 0], zlug, results_clay.get('Hinge location'))

    # # === SAND ===
    # y_sand, z_sand, results_sand = getCapacityDriven(profile_sand, 'sand', L, D, zlug, H=2.5e6, V=2.0e6, plot=True)

    # plot_pile(profile_sand, 'sand',  y_sand, z_sand, D, L, profile_sand[0, 0], zlug, results_sand.get('Hinge location'))

    # # === ROCK ===
    y_rock, z_rock, results_rock = getCapacityDriven(profile_rock, 'rock', L, D, zlug, H=3.5e6, V=3.0e6, plot=True)

    plot_pile(profile_rock, 'rock', y_rock, z_rock, D, L, profile_rock[0, 0], zlug, results_rock.get('Hinge location'))



