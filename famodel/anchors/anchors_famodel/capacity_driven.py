
import numpy as np
import matplotlib.pyplot as plt
from .support_soils import clay_profile, sand_profile, rock_profile
from .support_solvers import fd_solver
from .support_pycurves import py_Matlock, py_API, py_Reese
from .support_plots import plot_pile, plot_pycurve

def getCapacityDriven(profile_map, location_name, D, L, zlug, Ha, Va, plot=False):
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
    Ha : float
        Horizontal load applied at padeye (N)
    Va : float
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
    
    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']

    n = 50; loc = 2              # Number of nodes (-)
    tol = 1e-16; max_iter = 100  # Iteration parameters (-)
    nhuc = 1; nhu = 0.3          # Resistance factor (-)
    delta_r = 0.08               # Mean roughness height (m)
    
    t = (6.35 + D*20)/1e3        # Pile wall thickness (m), API RP2A-WSD
    E = 200e9                    # Elastic modulus of pile material (Pa)
    fy = 350e6                   # Steel's yield strength (Pa)
    rhows = 66.90e3              # Submerged steel specific weight (N/m3)
    rhow = 10e3                  # Water specific weight (N/m3) 

    I = (np.pi/64.0)*(D**4 - (D - 2*t)**4)
    EI = E*I
    h = L/n                      # Element size
    N = (n + 1) + 4              # (n+1) Real + 4 Imaginary nodes

    def PileSurface(Len, Dia):
        return np.pi*Dia*Len

    def PileWeight(Len, Dia, tw, rho):
        return ((np.pi/4)*(Dia**2 - (Dia - 2*tw)**2)*Len)*rho

    def SoilWeight(Len, Dia, tw, gamma_soil):
        return (np.pi/4)*(Dia - 2*tw)**2*Len*gamma_soil

    y = np.ones(N)*(0.01*D)
    z = np.zeros(N)
    py_funs  = []; PileShaft = []
    k_secant = np.zeros(N)
    DQ = []; pycurve_data = []
    
    z0 = min(layer['top'] for layer in layers)

    for i in [0, 1]:
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0

    for i in range(2, n+3):
        z[i] = (i - 2)*h
        z_depth = z[i]

        # Match the soil layer
        matched_layer = next((layer for layer in layers if layer['top'] <= z_depth <= layer['bottom']), None)
        if matched_layer is None or z_depth < matched_layer['top']:
            py_funs.append(lambda y_val: np.zeros_like(y_val))
            k_secant[i] = 0.0
            PileShaft.append(0.0)
            continue

        soil_type = matched_layer['soil_type']

        if soil_type == 'clay':
            profile = [[matched_layer['top'],    matched_layer['Su_top'], matched_layer['gamma_top']],
                       [matched_layer['bottom'], matched_layer['Su_bot'], matched_layer['gamma_bot']]]
            z0_local, f_Su, f_sigma_v_eff, f_gamma, f_alpha = clay_profile(profile)
            if z_depth < z0_local:
                py_funs.append(lambda y_val: np.zeros_like(y_val))
                k_secant[i] = 0.0
                PileShaft.append(0.0)
                continue
            Su = f_Su(z_depth)
            sigma_v_eff = f_sigma_v_eff(z_depth)
            gamma = f_gamma(z_depth)
            alpha = f_alpha(z_depth)
            py_f, (y_vals, p_vals) = py_Matlock(z_depth, D, gamma, Su, sigma_v_eff, z0=z0_local, return_curve=True)
            py_funs.append(py_f)
            pycurve_data.append((y_vals, p_vals, z_depth, 'clay'))
            Vo = np.pi*D*alpha*Su*z_depth**2
            PileShaft.append(Vo)
            k_val = py_funs[i](y[i])
            k_secant[i] = k_val/y[i] if y[i] != 0 else 0.0

        elif soil_type == 'sand':
            profile = [[matched_layer['top'],    matched_layer['phi_top'], matched_layer['gamma_top'], matched_layer['Dr_top']],
                       [matched_layer['bottom'], matched_layer['phi_bot'], matched_layer['gamma_bot'], matched_layer['Dr_bot']]]
            z0_local, f_phi, f_sigma_v_eff, f_gamma, f_Dr, f_delta = sand_profile(profile)
            if z_depth < z0_local:
                py_funs.append(lambda y_val: np.zeros_like(y_val))
                k_secant[i] = 0.0
                PileShaft.append(0.0)
                continue
            phi = f_phi(z_depth)
            sigma_v_eff = f_sigma_v_eff(z_depth)
            Dr = f_Dr(z_depth)
            delta = f_delta(z_depth)
            py_f, (y_vals, p_vals) = py_API(z_depth, D, phi, sigma_v_eff, Dr, z0=z0_local, return_curve=True)
            py_funs.append(py_f)
            pycurve_data.append((y_vals, p_vals, z_depth, 'sand'))
            fs = delta * sigma_v_eff
            Vo = np.pi*D*fs*z_depth
            PileShaft.append(Vo)
            k_val = py_funs[i](y[i])
            k_secant[i] = k_val/y[i] if y[i] != 0 else 0.0

        elif soil_type in ['rock', 'weak_rock']:
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
            py_f, (y_vals, p_vals) = py_Reese(z_depth, D, UCS, Em, z0=z0_local, return_curve=True)
            py_funs.append(py_f)
            pycurve_data.append((y_vals, p_vals, z_depth, 'rock'))
            SCR = nhuc*Em/(UCS*(1 + nhu))*delta_r/D
            alpha = 0.36*SCR - 0.0005
            fs = alpha*UCS
            Dq = np.pi*D*fs*z_depth
            DQ.append(Dq)
            k_val = py_funs[i](y[i])
            k_secant[i] = k_val/y[i] if y[i] != 0 else 0.0
            
        else:
            raise ValueError(f"Unsupported soil type: {matched_layer['soil_type']}")

    for i in [n+3, n+4]:
        z[i] = (i - 2)*h
        py_funs.append(0)
        k_secant[i] = 0.0

    gammas = [layer['gamma_top'] for layer in layers if layer['soil_type'] in ['clay', 'sand']]
    gamma_avg = np.mean(gammas) if gammas else 0
    Wp = PileWeight(L, D, t, rhows)
    Wsoil = SoilWeight(L, D, t, gamma_avg)
    Wshaft = PileShaft[-1] if PileShaft else 0.0
    last_layer_type = layers[-1]['soil_type']
    Wtip = DQ[-1] if DQ and last_layer_type in ['rock', 'weak_rock'] else 0.0

    Vmax = Wp + Wsoil + Wshaft + Wtip
                
    for j in range(max_iter):
        y_old = y.copy()
        y, Mi, Mp, hinge_formed, hinge_location = fd_solver(n, N, h, D, t, fy, EI, Ha, Va, zlug, z0, k_secant)
        
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
    y_pile = y[2:-2]
    z_pile = z[2:-2]
    ymax_index = np.argmax(np.abs(y_pile))

    resultsDriven = {
        'Horizontal max.': abs(Mi)/abs(zlug) if zlug != 0 else 1e-6,
        'Vertical max.': Vmax,
        'Lateral displacement': y_pile[ymax_index],
        'Rotational displacement': np.rad2deg(abs(y_pile[ymax_index - 1] - y_pile[ymax_index])/h),
        'Bending moment': abs(Mi),
        'Plastic moment': Mp,
        'Plastic hinge': hinge_formed,
        'Hinge location': hinge_location,         
        'Unity check (vertical)': Va/Vmax if Vmax != 0 else np.inf,
        'Unity check (horizontal)': Ha/(abs(Mi)/abs(zlug)) if zlug != 0 else np.inf,
        'Weight pile': PileWeight(L, D, t, rhows + rhow)}
    
    print(f"Max lateral displacement: {y_pile[ymax_index]:.6f} m at z = {z_pile[ymax_index]:.2f} m")
    print(f"Deflected tip: {y_pile[-1]:.6f} m at z = {z_pile[-1]:.2f} m")

    return layers, y[2:-2], z[2:-2], resultsDriven

if __name__ == '__main__':
    
    profile_map = [
        {
            'name': 'CPT_1',
            'x': 498234, 'y': 5725141,
            'layers': [
                {
                    'top': 1.0, 'bottom': 6.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.0, 'gamma_bot': 8.0,
                    'Su_top': 60, 'Su_bot': 200},
                # {
                #     'top': 6.0, 'bottom': 15.0,
                #     'soil_type': 'clay',
                #     'gamma_top': 8.0, 'gamma_bot': 8.0,
                #     'Su_top': 200, 'Su_bot': 400},
                {
                    'top': 6.0, 'bottom': 15.0,
                    'soil_type': 'sand',
                    'gamma_top': 8.0, 'gamma_bot': 8.0,
                    'phi_top': 32, 'phi_bot': 38,
                    'Dr_top': 70, 'Dr_bot': 75},
                {
                    'top': 15.0, 'bottom': 30.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.0, 'gamma_bot': 9.0,
                    'Su_top': 200, 'Su_bot': 400}]
        }
    ]

    D = 2.5           # Diameter (m)
    L = 25.0          # Length (m)
    zlug = 3          # Padeye depth (m)
    Ha = 5.0e5        # Horizontal load (N)
    Va = 1.5e5        # Vertical load (N)
    
    layers, y, z, resultsDriven = getCapacityDriven(profile_map, 'CPT_1', D, L, zlug, Ha, Va, plot=True)
    for key, val in resultsDriven.items():
        if isinstance(val, float):
            print(f"{key}: {val:.3f}")
        else:
            print(f"{key}: {val}")
    
    plot_pile(layers, y, z, D, L, z0=layers[0]['top'], zlug=zlug, hinge_location=None)




