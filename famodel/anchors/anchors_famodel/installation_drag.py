
import numpy as np
import matplotlib.pyplot as plt
from support_soils import clay_profile


def getInstallationDrag(profile_map, location_name, params, plot=True):
    '''Marching scheme for plate anchor trajectory and capacity.

    Parameters
    ----------
    profile : array
        Soil profile as a 2D array: (z, parameters)
            Clay soil profile (z (m), Su (kPa), gamma (kN/m³))
    location_name : str
        Name of the location in profile_map (e.g. 'CPT_1')
    params : dict
        Anchor and installation parameters (Af, Nc, En, b, mu, Ne,
        x0, z0, alpha, drag_install, drag_load, ds, q0i_deg, q0L_deg, fluke_shank_deg).
    n_steps : int
        Number of marching steps.

    Returns
    -------
    s : array
        Drag distance along seabed (m).
    depth : array
        Depth (m, positive down).
    Ta : array
        Mudline tension (N).
    qf_deg : array
        Fluke angle (deg).
    '''

    # --------------------------------------
    # 1) Retrieve soil layer for this location
    # --------------------------------------
    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']

    # Build combined clay profile: (z, gamma, Su) at all layer boundaries
    profile = []
    for layer in layers:
        z_top = layer['top']
        z_bot = layer['bottom']
        profile.append([z_top, layer['gamma_top'], layer['Su_top']])
        profile.append([z_bot, layer['gamma_bot'], layer['Su_bot']])

    # Call clay_profile exactly as in capacity_suction
    z_ref, f_gamma, f_Su, f_sigma_v_eff, f_alpha = clay_profile(profile)

    depth_min = profile[0][0]
    depth_max = profile[-1][0]

    # Helper: Su(z) and local gradient dSu/dz
    def su_and_k_depth(z):
        z = np.clip(z, depth_min, depth_max)
        Su = f_Su(z)
        k = (f_Su(np.clip(z + 0.1, depth_min, depth_max)) - f_Su(np.clip(z - 0.1, depth_min, depth_max)))/(0.2)
        return Su, k

    # --------------------------------------
    # 2) Read non-soil parameters from params
    # --------------------------------------
    Af = params['Af']
    Nc = params['Nc']
    En = params['En']
    b = params['b']
    mu = params['mu']          # not used explicitly in this sheet
    Ne = params['Ne']

    # x0, z0 are geometric / installation parameters (not soil)
    x0 = params['x0']          # initial drag distance along seabed (m)
    z0 = params['z0']          # initial embedment depth (m, positive down)

    alpha = params['alpha']
    drag_install = params['drag_install']
    drag_load = params['drag_load']
    ds = params['ds']
    q0i = np.radians(params['q0i_deg'])
    q0L = np.radians(params['q0L_deg'])
    fluke_shank = np.radians(params['fluke_shank_deg'])

    # Constants from Excel
    C = b*Nc*En/Af
    R = 0.01

    # --------------------------------------
    # 3) Allocate arrays
    # --------------------------------------
    n_steps=2000
    s = np.zeros(n_steps)       # drag distance
    z = np.zeros(n_steps)       # embedment (negative down)
    qa = np.zeros(n_steps)      # line angle at anchor
    qf = np.zeros(n_steps)      # fluke angle
    q0 = np.zeros(n_steps)      # mudline angle
    Su = np.zeros(n_steps)
    Su_avg = np.zeros(n_steps)
    Q = np.zeros(n_steps)       # operative bearing factor
    S = np.zeros(n_steps)       # dqadz
    T = np.zeros(n_steps)       # dx
    U = np.zeros(n_steps)       # dz
    W = np.zeros(n_steps)       # dqa
    V = np.zeros(n_steps)       # dqf
    Ta = np.zeros(n_steps)      # mudline tension (MN)

    # --------------------------------------
    # 4) Initial state (row 4 in Excel)
    # --------------------------------------
    s[0] = x0
    z[0] = -z0                  # negative down, consistent with original code
    q0[0] = q0i

    # Su at tip and average Su over initial embedment
    depth0 = -z[0]              # positive depth
    Su[0], k_eff0 = su_and_k_depth(depth0)
    Su_avg[0] = f_Su(0.5*depth0)
    Q[0] = Ne

    qa[0] = np.sqrt(2*b*Nc*En*abs(z[0])*Su_avg[0]/(Q[0]*Af*Su[0]) + q0[0]*q0[0])
    qf[0] = fluke_shank - qa[0]

    S[0] = (C/Q[0] - (qa[0]*qa[0] - q0[0]*q0[0])*k_eff0/(2.0*Su[0]))/qa[0]
    T[0] = ds*(np.cos(qf[0]) + R*np.sin(qf[0]))
    U[0] = -ds*(np.sin(qf[0]) - R*np.cos(qf[0]))
    W[0] = -U[0]*S[0]
    V[0] = W[0]

    # Initial mudline tension (converted to MN)
    Ta[0] = (Af*Q[0]*Su[0])*np.exp(alpha*qa[0])

    # --------------------------------------
    # 5) Marching loop
    # --------------------------------------
    for i in range(n_steps - 1):
        # 1) Update drag distance and embedment
        if s[i] < drag_install + drag_load:
            s_next = s[i] + T[i]
            z_next = z[i] + U[i]
        else:
            s_next = s[i]
            z_next = z[i]

        # 2) Mudline angle based on updated drag distance
        if s_next < drag_install:
            q0_next = q0i
        else:
            q0_next = q0L

        # 3) Fluke and line angles (B21 = 0 case → dqf = previous dqa)
        V_next = W[i]
        qf_next = qf[i] - V_next
        qa_next = qa[i] + W[i]

        # 4) Soil strength at new depth (still linear for now)
        depth_next = -z_next
        Su_next, k_eff_next = su_and_k_depth(depth_next)
        Su_avg_next = f_Su(0.5*depth_next)
        Q_next = Ne

        # 5) dqadz at next state
        S_next = (C/Q_next - (qa_next*qa_next - q0_next*q0_next)*k_eff_next/(2.0*Su_next))/qa_next

        # 6) Geometry increments for next step
        T_next = ds*(np.cos(qf_next) + R*np.sin(qf_next))
        U_next = -ds*(np.sin(qf_next) - R*np.cos(qf_next))
        W_next = -U_next*S_next

        # 7) Store next state
        s[i + 1] = s_next
        z[i + 1] = z_next
        q0[i + 1] = q0_next
        qf[i + 1] = qf_next
        qa[i + 1] = qa_next
        Su[i + 1] = Su_next
        Su_avg[i + 1] = Su_avg_next
        Q[i + 1] = Q_next
        S[i + 1] = S_next
        T[i + 1] = T_next
        U[i + 1] = U_next
        W[i + 1] = W_next
        V[i + 1] = V_next

        # 8) Mudline tension update (MN)
        if s_next < drag_install + drag_load:
            Ta[i + 1] = (Af*Q_next*Su_next)*np.exp(qa_next*alpha)
        else:
            Ta[i + 1] = Ta[i]

    depth = -z
    qf_deg = np.degrees(qf)
    
    if plot:
        # 1) Trajectory: drag distance vs depth
        plt.figure(figsize=(8, 6), dpi=100)
        plt.plot(s, depth, color='red', linestyle='-')
        plt.xlabel('Drag distance (m)')
        plt.ylabel('Depth (m)')
        plt.title('Anchor trajectory')
        plt.legend(['Drag Embedded Anchor (DEA) _ Clay Soil'], loc='upper left')
        plt.gca().invert_yaxis()
        plt.grid(True, which='major', linestyle='--')

        # 2) Fluke angle: drag distance vs fluke angle
        plt.figure(figsize=(8, 6), dpi=100)
        plt.plot(s, qf_deg, color='blue', linestyle='-')
        plt.xlabel('Drag distance (m)')
        plt.ylabel('Fluke angle (deg)')
        plt.legend(['Drag Embedded Anchor (DEA) _ Clay Soil'], loc='upper left')
        plt.title('Fluke angle vs drag distance')
        plt.grid(True, which='major', linestyle='--')

        # 3) Mudline tension: drag distance vs tension
        plt.figure(figsize=(8, 6), dpi=100)
        plt.plot(s, Ta, color='green', linestyle='-')
        plt.xlabel('Drag distance (m)')
        plt.ylabel('Mudline tension (N)')
        plt.legend(['Drag Embedded Anchor (DEA) _ Clay Soil'], loc='upper left')
        plt.title('Mudline tension vs drag distance')
        plt.grid(True, which='major', linestyle='--')
        
    resultsDrag = {}
    resultsDrag['Capacity'] = Ta[-1]
    resultsDrag['embedment_depth'] = depth[-1]
    resultsDrag['drag_distance'] = s[-1]
    resultsDrag['plate_angle'] = qf_deg[-1]
    
    return resultsDrag

if __name__ == '__main__':
    '''
    Testing the function with the same parameters
    as the Excel `trajectory-capacity` sheet.
    '''

    profile_map = [
        {
            'name': 'CPT_1',
            'x': 498234, 'y': 5725141,
            'layers': [
                {
                    'top': 0.0, 'bottom': 50.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.5, 'gamma_bot': 8.5,
                    'Su_top': 2.0, 'Su_bot': 102.0
                },
                # {
                #     'top': 10.0, 'bottom': 50.0,
                #     'soil_type': 'clay',
                #     'gamma_top': 8.5, 'gamma_bot': 8.5,
                #     'Su_top': 22.0, 'Su_bot': 102.0
                # }
            ]
        }
    ]

    params = {
        'Af': 19.64,
        'Lf': 4.43,          # not used directly, kept for completeness
        'b': 0.111,
        'Nc': 10.0,
        'En': 2.5,
        'mu': 0.5,
        'fluke_shank_deg': 42.0,
        'Ne': 7.8,
        'x0': 0.0,
        'z0': 1.0,
        'drag_install': 500.0,
        'q0i_deg': 0.0,
        'q0L_deg': 30.0,
        'drag_load': 0.0,
        'ds': 0.052,
        'alpha': 0.5,
    }

    resultsDrag = getInstallationDrag(profile_map, 'CPT_1', params, plot=True)
