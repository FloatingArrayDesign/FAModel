
import numpy as np
import matplotlib.pyplot as plt
from .support_soils import rock_profile
from .support_solvers import fd_solver
from .support_pycurves import py_Lovera
from .support_plots import plot_pile, plot_pycurve

def getCapacityDrilled(profile_map, location_name, L, D, zlug, Ha, Va, plot=False, display=0):
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

    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']
      
    n = 50; loc = 2              # Number of nodes (-)
    tol = 1e-16; max_iter = 50   # Iteration parameters (-)
    nhuc = 1; nhu = 0.3          # Resistance factor (-)
    delta_r = 0.08               # Mean roughness height (m)
    
    t = (6.35 + D*20)/1e3        # Pile wall thickness (m), API RP2A-WSD
    E = 200e9                    # Elastic modulus of pile material (Pa)
    fy = 350e6                   # Steel's yield strength (Pa)
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
    DQ = []; pycurve_data = []

    z0 = min(layer['top'] for layer in layers)
    
    for i in [0, 1]:              # Top two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0
   
    for i in range(2, n+3):       # Real nodes
        z[i] = (i - 2)*h
        z_depth = z[i]
        
        matched_layer = next((layer for layer in layers if layer['top'] <= z_depth <= layer['bottom']), None)
        if matched_layer is None or z_depth < matched_layer['top']:
            py_funs.append(lambda y_val: np.zeros_like(y_val))
            k_secant[i] = 0.0
            DQ.append(0.0)
            continue
        
        profile = [[matched_layer['top'],    matched_layer['UCS_top'], matched_layer['Em_top']],
                   [matched_layer['bottom'], matched_layer['UCS_bot'], matched_layer['Em_bot']]]
        z0_local, f_UCS, f_Em = rock_profile(profile)
        
        if z_depth < z0_local:
            py_funs.append(lambda y_val: np.zeros_like(y_val))
            k_secant[i] = 0.0
            DQ.append(0.0)
            continue

        UCS = f_UCS(z_depth)
        Em = f_Em(z_depth)
        py_f, (y_vals, p_vals) = py_Lovera(z_depth, D, UCS, Em, zlug, z0, return_curve=True)
        py_funs.append(py_f)
        pycurve_data.append((y_vals, p_vals, z_depth, 'rock'))
        # print(f"z_depth = {z_depth:.2f} m, UCS = {f_UCS(z_depth):.2e} Pa, Em = {f_Em(z_depth):.2e} Pa")
        
        SCR = nhuc*Em/(UCS*(1 + nhu))*delta_r/D
        alpha = 0.36*SCR - 0.0005
        fs = alpha*UCS
        Dq = np.pi*D*fs*z_depth
        DQ.append(Dq)
        k_val = py_funs[i](y[i])
        k_secant[i] = k_val/y[i] if y[i] != 0 else 0.0

    for i in [n+3, n+4]:         # Bottom two imaginary nodes
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0
        
    Wp = PileWeight(L, D, t, rhows + rhow)
    Wtip = DQ[-1] if DQ else 0.0
    Vmax = Wp + Wtip    
        
    for j in range(max_iter):
        y_old = y.copy()
        y, *_ = fd_solver(n, N, h, D, t, fy, EI, Ha, Va, zlug, z0, k_secant)
        
        # Update stiffness
        for i in range(2, n+3):
            if callable(py_funs[i]):
                k_secant[i] = py_funs[i](y[i])/y[i] if y[i] != 0 else 0.0
    
        # Check convergence
        if np.linalg.norm(y - y_old, ord=2) < tol:
            if display > 0: print(f'[Converged in {j+1} iterations]')
            break
    else:
        if display > 0: print('[Warning: Solver did not converge]')


    if plot:
        plot_pycurve(pycurve_data)

        fig, ax = plt.subplots(figsize=(3, 5))
        y0 = np.zeros_like(z[2:-2])
        ax.plot(y0, z[2:-2], 'k', label='Original pile axis')
        ax.plot(y[2:-2], z[2:-2], 'r', label='Deflected shape')
        ax.plot(0, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')
        ax.axhline(z0, color='blue', linestyle='--', label=f'Mudline (z0 = {z0:.2f} m)')
        ax.set_xlabel('Lateral displacement (m)')
        ax.set_ylabel('Depth (m)')
        ax.set_xlim([-0.1*D, 0.1*D])
        ax.set_ylim([L + 5, -2])
        ax.grid(ls='--')
        ax.legend() 
            
    # Relevant index of nodes
    zlug_index = int(zlug/h)
    if display > 0: print(zlug_index)
    ymax_index = np.argmax(y)
    if display > 0: print(ymax_index)    
    
    resultsDandG = {
        'Vertical max.': Vmax,
        'Lateral displacement': y[ymax_index],
        'Rotational displacement': np.rad2deg(abs(y[ymax_index - 1] - y[ymax_index])/h),
        'Bending moment': None,
        'Plastic moment': None,
        'Plastic hinge': None,
        'Hinge location': None,      
        'Unity check (vertical)': Va/Vmax if Vmax != 0 else np.inf,
        'Unity check (horizontal)': 0.0,  # Placeholder; no Mp or Mi in current model
        'Weight pile': PileWeight(L, D, t, rhows + rhow),
        'p-y model': 'Lovera (2023)'}
    
    return layers, y[2:-2], z[2:-2], resultsDandG
   
if __name__ == '__main__':

    profile_map = [
        {
            'name': 'CPT_rock_1',
            'x': 502000,
            'y': 5725000,
            'layers': [
                {
                    'top': 2.0, 'bottom': 5.0,
                    'soil_type': 'rock',
                    'UCS_top': 8.0, 'UCS_bot': 8.0,    # MPa
                    'Em_top': 100, 'Em_bot': 200       # MPa
                },
                {
                    'top': 5.0, 'bottom': 9.0,
                    'soil_type': 'rock',
                    'UCS_top': 10.0, 'UCS_bot': 10.0,  # MPa
                    'Em_top': 200, 'Em_bot': 300       # MPa
                },
                {
                    'top': 9.0, 'bottom': 30.0,
                    'soil_type': 'rock',
                    'UCS_top': 20.0, 'UCS_bot': 20.0,  # MPa
                    'Em_top': 300, 'Em_bot': 400       # MPa
                }
            ]
        }
    ]

    D = 3.0           # Diameter (m)
    L = 10.0          # Length (m)
    zlug = 1          # Padeye elevation (m)
    Ha = 5.0e6        # Horizontal load (N)
    Va = 3.0e5        # Vertical load (N)

    layers, y, z, results = getCapacityDrilled(profile_map, 'CPT_rock_1', L, D, zlug, Ha, Va, plot=True, display=0)

    print('\n--- Results for Drilled Pile in Layered Rock ---')
    for key, val in results.items():
        print(f'{key}: {val:.3f}' if isinstance(val, float) else f'{key}: {val}')
        
    plot_pile(layers, y, z, D, L, layers[0]['top'], zlug)    

    

    
    
