import numpy as np
import matplotlib.pyplot as plt
from capacity_soils_map import clay_profile, sand_profile, rock_profile

def getInstallationDriven(profile_map, location_name, D, L, hammer, J_shaft, J_toe, plot=True, refusal_threshold=0.002, refusal_count=10):
    dz = 0.24
    z0 = 1.0
    max_depth = L
    N = int((max_depth - z0) / dz)

    t_wall = (6.35 + D * 20) / 1e3
    D_inner = D - 2 * t_wall
    A = np.pi / 4 * (D**2 - D_inner**2)
    E = 2.1e11
    rhos = 7850
    g = 9.81

    m = rhos*A*dz
    k = E*A/dz
    dt = N/np.sqrt(E/rhos)

    m_r = hammer['m_r']
    h = hammer['h']
    eta = hammer['efficiency']
    E_hammer = eta*m_r*g*h
    v0 = np.sqrt(2 * E_hammer / m_r)

    soil = profile_map[location_name]
    layers = soil['layers']

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

    penetration = 0.0
    total_energy = 0.0
    toe_displacement = []
    blow_counts = []
    consecutive_small_blows = 0

    blow = 0
    while penetration < (max_depth - z0):
        u = np.zeros(N + 1)
        v = np.zeros(N + 1)
        a = np.zeros(N + 1)

        # Apply initial velocity to first node to help energy propagate
        v[0] = v0
        T = 0.5
        nsteps = int(T/dt)

        for step in range(nsteps):
            F_internal = np.zeros(N + 1)
            for i in range(1, N):
                F_internal[i] = k*(u[i - 1] - 2*u[i] + u[i + 1])

            F_internal[N] = k*(u[N - 1] - u[N])

            for i in range(N + 1):
                z = dz*i + z0
                
                if i == N:
                    F_internal[i] -= compute_Rdynamic(v[i], z, J_toe)
                    F_internal[i] += compute_Rstatic(z)
                else:
                    F_internal[i] -= compute_Rdynamic(v[i], z, J_shaft)

                if step < 10:
                    print(f"  Step {step}: z = {z:.2f}, Rd = {compute_Rdynamic(v[i], z, J_shaft):.2f}, v = {v[i]:.4f}, u = {u[i]:.4f}")

            a = np.nan_to_num(F_internal / m, nan=0.0, posinf=0.0, neginf=0.0)
            v += a * dt
            u += v * dt

            if np.any(np.abs(u) > 5):
                print("Displacement blew up. Stopping.")
                break

        delta_z = u[-1]
        penetration += delta_z

        toe_displacement.append(penetration)
        blow_counts.append(blow + 1)
        total_energy += E_hammer

        print(f"Blow {blow + 1}: Δz = {delta_z:.5f} m, Total Penetration = {penetration:.3f} m")

        if abs(delta_z) < refusal_threshold:
            consecutive_small_blows += 1
        else:
            consecutive_small_blows = 0

        if consecutive_small_blows >= refusal_count:
            print("Refusal criteria met: 10 consecutive blows with <2 mm displacement")
            break

        blow += 1

    if plot:
        plt.figure(figsize=(8, 4))
        plt.plot(blow_counts, toe_displacement, marker='o')
        plt.xlabel('Blow Count')
        plt.ylabel('Cumulative Toe Displacement (m)')
        plt.title('Toe Displacement vs Blow Count')
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return {
        'blow_counts': blow_counts,
        'toe_displacement': toe_displacement,
        'total_energy': total_energy,
        'final_depth': penetration,
        'total_counts': len(blow_counts)
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
    hammer = {'m_r': 85000, 'h': 5.5, 'efficiency': 0.85}
    J_shaft = 0.05
    J_toe = 0.05

    results = getInstallationDriven(profile_map, 'CPT_1', D, L, hammer, J_shaft, J_toe, plot=True)
    for key, val in results.items():
        print(f"{key}: {val}")
