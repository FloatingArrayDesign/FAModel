import numpy as np
import matplotlib.pyplot as plt
from capacity_soils_map import clay_profile, sand_profile, rock_profile

def compute_Rstatic(z):
    for layer in layers:
        if layer['top'] <= z <= layer['bottom']:
            if layer['soil_type'] == 'clay':
                profile = [[layer['top'],    layer['gamma_top'], layer['Su_top']],
                           [layer['bottom'], layer['gamma_bot'], layer['Su_bot']]]
                _, _, _, _, f_alpha = clay_profile(profile)
                return f_alpha(z)*np.pi*D*dz
            elif layer['soil_type'] == 'sand':
                profile = [[layer['top'],    layer['gamma_top'], layer['phi_top'], layer['Dr_top']],
                           [layer['bottom'], layer['gamma_bot'], layer['phi_bot'], layer['Dr_bot']]]
                _, _, _, _, _, f_delta = sand_profile(profile)
                return f_delta(z)*np.pi*D*dz
            elif layer['soil_type'] == 'rock':
                profile = [[layer['top'],    layer['UCS_top'], layer['Em_top']],
                           [layer['bottom'], layer['UCS_bot'], layer['Em_bot']]]
                _, f_UCS, _ = rock_profile(profile)
                return f_UCS(z)*np.pi*D*dz
    return 0.0

def compute_Rdynamic(v_local, z, J):
    for layer in layers:
        if layer['top'] <= z <= layer['bottom']:
            if layer['soil_type'] == 'clay':
                profile = [[layer['top'],    layer['gamma_top'], layer['Su_top']],
                           [layer['bottom'], layer['gamma_bot'], layer['Su_bot']]]
                _, _, _, _, f_alpha = clay_profile(profile)
                return J*f_alpha(z)*v_local*np.pi*D
            elif layer['soil_type'] == 'sand':
                profile = [[layer['top'],    layer['gamma_top'], layer['phi_top'], layer['Dr_top']],
                           [layer['bottom'], layer['gamma_bot'], layer['phi_bot'], layer['Dr_bot']]]
                _, _, _, _, _, f_delta = sand_profile(profile)
                return J*f_delta(z)*v_local*np.pi*D
            elif layer['soil_type'] == 'rock':
                profile = [[layer['top'],    layer['UCS_top'], layer['Em_top']],
                           [layer['bottom'], layer['UCS_bot'], layer['Em_bot']]]
                _, f_UCS, _ = rock_profile(profile)
                return J*f_UCS(z)*v_local*np.pi*D
    return 0.0

def getInstallationDriven(profile_map, location_name, D_input, L, hammer, J_shaft, J_toe, plot):
    global D, dz, layers
    D = D_input
    soil = profile_map[location_name]
    layers = soil['layers']
    z0 = layers[0]['top']

    N = 100
    z = np.linspace(z0, L, N)
    dz = z[1] - z[0]

    dt = 0.001
    t_max = 2.0
    time = np.arange(0, t_max, dt)

    u = np.zeros((len(time), N))
    v = np.zeros((len(time), N))
    F = np.zeros((len(time), N))

    rho_steel = 7850
    t = 0.05
    A = np.pi * ((D / 2)**2 - ((D / 2) - t)**2)
    m_total = A * L * rho_steel
    M_prime = m_total / N

    m_r = hammer['m_r']
    h = hammer['h']
    eff = hammer['efficiency']
    E = m_r * 9.81 * h * eff

    blow_count = 0
    penetration = 0.0
    z_vals = []
    blow_vals = []
    refusal_counter = 0
    refusal_limit = 10
    min_set = 0.002

    while penetration < L and refusal_counter < refusal_limit:
        blow_count += 1
        v0 = np.sqrt(2 * E / m_r)
        v[0, 0] = v0

        for n in range(1, len(time)):
            for i in range(1, N - 1):
                u_xx = (u[n - 1, i + 1] - 2 * u[n - 1, i] + u[n - 1, i - 1]) / dz**2
                F[n, i] = M_prime * u_xx

            v[n] = v[n - 1] + (F[n] / M_prime) * dt
            u[n] = u[n - 1] + v[n] * dt

        delta_z = max(u[-1, :])
        tip_depth = penetration + delta_z

        Rs = compute_Rstatic(tip_depth)
        Rt = 9.0 * compute_Rstatic(tip_depth)
        R_total = Rs + Rt + compute_Rdynamic(v[-1, -1], tip_depth, J_shaft)

        penetration += delta_z
        z_vals.append(penetration)
        blow_vals.append(blow_count)

        if delta_z < min_set:
            refusal_counter += 1
        else:
            refusal_counter = 0

    if plot:
        plt.figure()
        plt.plot(z_vals, blow_vals, label='Blows vs Penetration')
        plt.xlabel('Penetration Depth [m]')
        plt.ylabel('Blow Count')
        plt.title('Driveability Simulation')
        plt.grid(True)
        plt.legend()
        plt.show()

    return {
        'z': z,
        'z_vals': z_vals,
        'blow_vals': blow_vals,
        'u': u,
        'v': v,
        'F': F,
        'dt': dt,
        'dz': dz
    }


if __name__ == '__main__':
    profile_map = {
        'CPT_1': {
            'type': 'clay',
            'x': 498234, 'y': 5725141,
            'layers': [
                {
                    'top': 1.0, 'bottom': 6.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.0, 'gamma_bot': 8.0,
                    'Su_top': 10, 'Su_bot': 20},
                {
                    'top': 6.0, 'bottom': 15.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.0, 'gamma_bot': 8.0,
                    'Su_top': 80, 'Su_bot': 100},
                {
                    'top': 15.0, 'bottom': 30.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.0, 'gamma_bot': 9.0,
                    'Su_top': 100, 'Su_bot': 200}
            ]
        }
    }

    D = 1.0
    L = 12.0
    hammer = {'m_r': 155000, 'h': 5.5, 'efficiency': 0.85}
    J_shaft = 0.05
    J_toe = 0.05

    results = getInstallationDriven(profile_map, 'CPT_1', D, L, hammer, J_shaft, J_toe, plot=True)
    for key, val in results.items():
        print(f"{key}: {val}")