
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import linalg
from .capacity_soils import rock_profile
from .capacity_pycurves import py_Lovera
from .capacity_plots import plot_pile

def getCapacityDandG(profile, soil_type, L, D, zlug, Ha, Va, plot=True):
    '''Models a laterally loaded pile using the p-y method. The solution for
    lateral displacements is obtained by solving the 4th order ODE, 
    EI*d4y/dz4 - V*d2y/dz2 + ky = 0 using the finite difference method.
    EI*d4y/dz4 - V*d2y/dz2 + K*z*dy/dz + ky = 0 using the finite difference method.

    Assumes that EI remains constant with respect to curvature i.e. pile
    material remains in the elastic region.

    Parameters
    ----------
    profile : array
        Rock profile as a 2D array: (z (m), UCS (MPa), Em (MPa))
    soil_type : string
        Select soil condition, 'rock'
    L : float 
        Pile length (m)
    D : float 
        Pile diameter (m)
    zlug : float
        Load eccentricity above the mudline or depth to mudline relative to the pile head (m)
    Ha : float       
        Horizontal load at pile lug elevation (N)
    Va : float          
        Vertical load at pile lug elevation (N)
    plot : bool
        Plot the p-y curve and the deflection pile condition if True

    Returns
    -------
    y : array
        Lateral displacement at each node (n+1 real + 4 imaginary)
    z : array
        Node location along pile (m)
    resultsDandG : dict
        Dictionary with lateral, rotational, vertical and pile weight results
    '''

    n = 50; loc = 2              # Number of nodes (-)
    tol = 1e-16; max_iter = 50   # Iteration parameters (-)
    nhuc = 1; nhu = 0.3          # Resistance factor (-)
    delta_r = 0.08               # Mean roughness height (m)
    
    t = (6.35 + D*20)/1e3        # Pile wall thickness (m), API RP2A-WSD
    E = 200e9                    # Elastic modulus of pile material (Pa)
    rhows = 66.90e3              # Submerged steel specific weight (N/m3)
    rhow = 10e3                  # Water specific weight (N/m3) 
    
    # Pile geometry
    I = (np.pi/64.0)*(D**4 - (D - 2*t)**4)
    EI = E*I
    h = L/n                      # Element size
    N = (n + 1) + 4              # (n+1) Real + 4 Imaginary nodes
    
    # Dry and wet mass of the pile    
    def PileWeight(Len, Dia, tw, rho):
        Wp = ((np.pi/4)*(Dia**2 - (Dia - 2*tw)**2)*Len)*rho
        return Wp 

    # Array for displacements at nodes, including imaginary nodes.
    y = np.ones(N)*(0.01*D)      # An initial value of 0.01D was arbitrarily chosen

    # Initialize and assemble array/list of p-y curves at each real node
    z = np.zeros(N)
    k_secant = np.zeros(N)
    py_funs = []
    DQ = []
    z0, f_UCS, f_Em = rock_profile(profile)
    

    for i in [0, 1]:              # Top two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0
   
    for i in range(2, n+3):       # Real nodes
        z[i] = (i - 2)*h       
        if z[i] < z0:
             # No p-y curve above mudline
             py_funs.append(lambda y_val: np.zeros_like(y_val))
             k_secant[i] = 0.0
             DQ.append(0.0)
        else:
             py_funs.append(py_Lovera(z[i], D, f_UCS, f_Em, zlug, z0, plot=True)) 
             # print(f"z = {z[i]:.2f} m, UCS = {f_UCS(z[i]):.2e} Pa, Em = {f_Em(z[i]):.2e} Pa")
             UCS = f_UCS(z[i])
             Em = f_Em(z[i])
             SCR = nhuc*Em/(UCS*(1 + nhu))*delta_r/D
             alpha = 0.36*SCR - 0.0005
             fs = alpha*UCS
             Dq = np.pi*D*fs*z[i]
             DQ.append(Dq)
             Vmax = PileWeight(L, D, t, rhows) + DQ[-1]
             
             k_secant[i] = py_funs[i](y[i])/y[i]

    for i in [n+3, n+4]:         # Bottom two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0
   
    for j in range(max_iter):
        y_old = y.copy()
        y = fd_solver(n, N, h, EI, Ha, Va, zlug, z0, k_secant)
        
        # Update stiffness
        for i in range(2, n+3):
            if callable(py_funs[i]):
                k_secant[i] = py_funs[i](y[i])/y[i] if y[i] != 0 else 0.0
    
        # Check convergence
        if np.linalg.norm(y - y_old, ord=2) < tol:
            print(f'[Converged in {j+1} iterations]')
            break
    else:
        print('[Warning: Solver did not converge]')

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
            
    # Relevant index of nodes
    zlug_index = int(zlug/h); print(zlug_index)
    ymax_index = int(np.max(y)); print(ymax_index)    
    
    resultsDandG = {
        'Weight pile': PileWeight(L, D, t, rhows + rhow),
        'Vertical max.': Vmax,
        'Lateral displacement': y[ymax_index],
        'Rotational displacement': np.rad2deg(abs(y[ymax_index - 1] - y[ymax_index])/h),
        'Unity check (vertical)': Va/Vmax if Vmax != 0 else np.inf,
        'Unity check (horizontal)': 0.0,  # Placeholder; no Mp or Mi in current model
        'Bending moment': None,
        'Plastic moment': None,
        'Plastic hinge': None,
        'Hinge location': None,
        'p-y model': 'Lovera (2023)',
        }
    
    return y[2:-2], z[2:-2], resultsDandG

def fd_solver(n, N, h, EI, Ha, Va, zlug, z0, k_secant):
    '''Solves the finite difference equations from 'py_analysis_1'. This function should be run iteratively for
    non-linear p-y curves by updating 'k_secant' using 'y'. A single iteration is sufficient if the p-y curves
    are linear.

    Parameters
    ----------
    n : int
        Number of elements (-)
    N : int
        Total number of nodes (-)
    h : float
        Element size (m)
    EI : float
        Flexural rigidity of the pile (Nm²)
    Ha : float
        Horizontal load at padeye (N)
    Va : float
        Vertical load at padeye (N)
    zlug : float
        Padeye depth from pile head (m)
    z0 : float
        Mudline elevation from pile head (m)
    k_secant : array
        Secant stiffness at each node

    Returns
    -------
    y : array
        Lateral displacement at each node (n+1 real + 4 imaginary)
    '''
    
    # Initialize and assemble matrix
    X = np.zeros((N, N))

    # (n+1) finite difference equations for (n+1) real nodes
    for i in range(0,n+1):
        X[i, i]   =  1.0
        X[i, i+1] = -4.0 + Va*h**2/EI
        X[i, i+2] =  6.0 - 2*Va*h**2/EI + k_secant[i+2]*h**4/EI
        X[i, i+3] = -4.0 + Va*h**2/EI
        X[i, i+4] =  1.0

    # Curvature at pile head
    X[n+1, 1] =  1.0
    X[n+1, 2] = -2.0
    X[n+1, 3] =  1.0

    # Shear at pile head
    X[n+2, 0] = -1.0
    X[n+2, 1] =  2.0 - Va*h**2/EI
    X[n+2, 2] =  0.0
    X[n+2, 3] = -2.0 + Va*h**2/EI
    X[n+2, 4] =  1.0

    # Curvature at pile tip
    X[n+3, -2] =  1.0
    X[n+3, -3] = -2.0
    X[n+3, -4] =  1.0

    # Shear at pile tip
    X[n+4, -1] =  1.0
    X[n+4, -2] = -2.0 + Va*h**2/EI
    X[n+4, -3] =  0.0
    X[n+4, -4] =  2.0 - Va*h**2/EI
    X[n+4, -5] = -1.0

    # Initialize vector q
    q = np.zeros(N)
    
    # Index of the node where the horizontal load is applied (padeye)
    zlug_index = int(zlug/h)
    q[zlug_index] = 2*Ha*h**3

    # Solve for displacement
    y = linalg.solve(EI*X, q)

    return y
   
if __name__ == '__main__':
    
    profile_rock = np.array([
        [ 2.0, 2, 200],
        [ 5.0, 2, 200],
        [ 9.0, 2, 200],
        [30.0, 2, 200]
    ])
    
    D = 3.0           # Diameter (m)
    L = 10.0          # Length (m)
    zlug = 1          # Padeye elevation (m)
    Ha = 8.0e6        # Horizontal load (N)
    Va = 3.0e6        # Vertical load (N)
    
    y, z, resultsDandG = getCapacityDandG(profile_rock, 'rock', L, D, zlug, Ha, Va, plot=True)
    for key, val in resultsDandG.items():
        if isinstance(val, float):
            print(f"{key}: {val:.3f}")
        else:
            print(f"{key}: {val}")
    
    plot_pile(profile_rock, 'rock', y, z, D, L, profile_rock[0, 0], zlug)
