
import numpy as np
import matplotlib.pyplot as plt
from .support_soils import clay_profile, sand_profile
from .support_plots import plot_load

def getTransferLoad(profile_map, Tm, thetam, zlug, line_type, d, w=None, plot=True):
    '''Calculate the transfer load from mudline to main padeye using a layered soil profile.

    Parameters
    ----------
    profile_map : list of dicts
        Soil profile in profile_map format
    Tm : float
        Mooring line load at mudlevel (N)
    thetam : float
        Mooring line angle at mudlevel (deg)
    zlug : float
        Embedment depth of the lug (m)
    line_type : str
        'chain' or 'wire'
    d : float
        Chain diameter (m)
    w : float
        Mooring line unit weight (N/m)
    plot : bool
        Show plot

    Returns
    -------
    dict
        Dictionary with transferred load components and depth.
    '''

    deltas = 0.2  # discretization step

    # Line mechanical properties
    if line_type == 'chain':
        Et, En = 10, 2.5
    elif line_type == 'wire':
        Et, En = np.pi, 1
    W = w*deltas

    # Soil layer access
    layers = profile_map[0]['layers']
    z0 = min(layer['top'] for layer in layers)
    Nc = 8.5

    # Initial values
    z0 = min(layer['top'] for layer in layers)
    T = Tm
    theta = np.deg2rad(thetam)
    drag = 0
    depth = z0 + 0.01

    # Tracing lists
    drag_values, depth_values = [], []

    while (zlug - depth) >= 0:
        matched_layer = next((layer for layer in layers if layer['top'] <= depth <= layer['bottom']), None)
        if matched_layer is None:
            break
        
        if matched_layer['soil_type'] == 'clay':
            matched_layer = next((layer for layer in layers if layer['soil_type'] == 'clay' and layer['top'] <= depth <= layer['bottom']), None)
            if matched_layer is None:
                break
            profile = [[matched_layer['top'],    matched_layer['gamma_top'], matched_layer['Su_top']],
                       [matched_layer['bottom'], matched_layer['gamma_bot'], matched_layer['Su_bot']]]
            z0_local, f_gamma, f_Su, f_sigma_v_eff, f_alpha = clay_profile(profile)

            Su = f_Su(depth)
            alpha = f_alpha(depth)
            d_theta = (En*d*Nc*Su - W*np.cos(theta))/T*deltas
            dT = (Et*d*alpha*Su + W*np.sin(theta))*deltas

        elif matched_layer['soil_type'] == 'sand':
            matched_layer = next((layer for layer in layers if layer['soil_type'] == 'sand' and layer['top'] <= depth <= layer['bottom']), None)
            if matched_layer is None:
                break
            
            profile = [[matched_layer['top'],    matched_layer['gamma_top'], matched_layer['phi_top'], matched_layer['Dr_top']],
                       [matched_layer['bottom'], matched_layer['gamma_bot'], matched_layer['phi_bot'], matched_layer['Dr_bot']]]
            z0_local, f_gamma, f_phi, f_Dr, f_sigma_v_eff, f_delta = sand_profile(profile)

            gamma_z = f_gamma(depth)
            delta_z = f_delta(depth)
            phi = f_phi(depth)
            Nq = np.exp(np.pi*np.tan(np.deg2rad(phi)))*(np.tan(np.deg2rad(45 + phi/2)))**2
            print(f'Nq = {Nq:.2f}, depth = {depth:.2f} m')
            d_theta = (En*d*Nq*gamma_z*depth - W*np.cos(theta))/T*deltas
            dT = (Et*d*gamma_z*depth*np.tan(np.deg2rad(delta_z)) + W*np.sin(theta))*deltas
            
        else:
            raise ValueError(f"Unsupported soil type: {matched_layer['soil_type']}")

        d_drag  = deltas*np.cos(theta)
        d_depth = deltas*np.sin(theta)

        theta += d_theta
        T -= dT
        drag += d_drag
        depth += d_depth

        if abs(Tm - T) > 0.75*Tm:
            raise Exception(f"Load transfer unrealistic: Tm = {Tm/1e6:.2f} MN vs T = {T/1e6:.2f} MN")
        if not (0 < np.rad2deg(theta) < 90):
            raise Exception(f"Load angle unrealistic: {np.rad2deg(theta):.2f} deg")

        drag_values.append(-drag); 
        depth_values.append(-depth); 

    Ta = T; thetaa = theta
    Hm = Tm*np.cos(np.deg2rad(thetam)); Vm = Tm*np.cos(np.deg2rad(thetam))
    Ha = Ta*np.cos(thetaa); Va = Ta*np.sin(thetaa)
    
    print(f'Input Tm = {Tm}, thetam = {thetam}, zlug = {zlug}')
    print(f'Output Hm = {Hm}, Vm = {Vm}')
    print(f'Output Ta = {Ta}, thetaa = {np.rad2deg(thetaa)}')
    print(f'Output Ha = {Ha}, Va = {Va}')

    resultsLoad = {
        'Tm': Tm, 'thetam': thetam,
        'Hm': Hm, 'Vm': Vm,
        'Ta': Ta, 'thetaa': np.rad2deg(thetaa),
        'Ha': Hm, 'Va': Vm,
        'length': deltas*len(drag_values),
        'drag_values': drag_values,
        'depth_values': depth_values}

    return layers, resultsLoad


if __name__ == '__main__':

    # profile_map = [
    #     {
    #         'name': 'CPT_1',
    #         'x': 498234, 'y': 5725141,
    #         'layers': [
    #             {
    #                 'top': 1.0, 'bottom': 2.0,
    #                 'soil_type': 'clay',
    #                 'gamma_top': 8.0, 'gamma_bot': 8.0,
    #                 'Su_top': 10, 'Su_bot': 25},
    #             {
    #                 'top': 2.0, 'bottom': 8.0,
    #                 'soil_type': 'clay',
    #                 'gamma_top': 8.0, 'gamma_bot': 8.0,
    #                 'Su_top': 25, 'Su_bot': 50},
    #             {
    #                 'top': 8.0, 'bottom': 16.0,
    #                 'soil_type': 'clay',
    #                 'gamma_top': 8.0, 'gamma_bot': 8.0,
    #                 'Su_top': 50, 'Su_bot': 100}
    #         ]
    #     }
    # ]
    profile_map = [
        {
            'name': 'CPT_1',
            'x': 498234, 'y': 5725141,
            'layers': [
                # {
                #     'top': 0.0, 'bottom': 5.0,
                #     'soil_type': 'sand',
                #     'gamma_top': 9.5, 'gamma_bot': 9.5,
                #     'phi_top': 28, 'phi_bot': 30,
                #     'Dr_top': 70, 'Dr_bot': 70},
                # {
                #     'top': 0.0, 'bottom': 5.0,
                #     'soil_type': 'clay',
                #     'gamma_top': 8.0, 'gamma_bot': 8.0,
                #     'Su_top': 25, 'Su_bot': 25},
                {
                    'top': 0.0, 'bottom': 3.0,
                    'soil_type': 'sand',
                    'gamma_top': 9.5, 'gamma_bot': 9.5,
                    'phi_top': 25, 'phi_bot': 30,
                    'Dr_top': 60, 'Dr_bot': 65},
                {
                    'top': 3.0, 'bottom': 15.0,
                    'soil_type': 'sand',
                    'gamma_top': 9.5, 'gamma_bot': 9.5,
                    'phi_top': 32, 'phi_bot': 35,
                    'Dr_top': 70, 'Dr_bot': 85}
            ]
        }
    ]

    Tm = 4978442          # Load at mudline (N)
    thetam = 15            # Angle at mudline (deg)
    zlug = 8.5              # Padeye depth (m)
    line_type = 'chain'
    d = 0.12              # Chain diameter (m)
    w = 2000              # Line weight (N/m)

    layers, resultsLoad = getTransferLoad(profile_map, Tm, thetam, zlug, line_type, d, w, plot=True)

    # print("\n--- Transfer Load Results ---")
    # for key, val in resultsLoad.items():
    #     if isinstance(val, float):
    #         print(f"{key}: {val:.3f}")
    #     elif isinstance(val, list):
    #         print(f"{key}:")
    #         for v in val:
    #             print(f"  {v:.3f}")
    #     else:
    #         print(f"{key}: {val}")

    plot_load(layers, resultsLoad['drag_values'], resultsLoad['depth_values'], 
              resultsLoad['Tm'], resultsLoad['thetam'], resultsLoad['Ta'], 
              resultsLoad['thetaa'], zlug=zlug)