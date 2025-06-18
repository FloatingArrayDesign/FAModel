
import numpy as np
import matplotlib.pyplot as plt
from support_soils import sand_profile

def getInstallationHelical(profile_map, location_name, D, L, d, plot=True):
    # Constants and geometry
    dH = 0.05                  # m
    psi_p = 16.5               # degrees
    delta_crit = 24            # degrees
    Fr = 0.01

    Dh = D                     # Helix diameter
    Dc = d                     # Core diameter
    ph = Dh/3                  # Pitch
    E = 210e9                  # Pa
    th = 0.10                  # m
    tc = 0.03                  # m
    f_y = 350e6                # Pa
    f_y_weld = 250e6           # Pa (typical weld yield strength)
    k = 1.04                   # Bending coefficient

    profile_entry = next((p for p in profile_map if p['name'] == location_name), None)
    if not profile_entry:
        raise ValueError(f"Location '{location_name}' not found in profile_map")

    layers = profile_entry['layers']

    # Assemble sand profile
    profile = []
    for layer in layers:
        if layer['soil_type'] == 'sand':
            profile.append([layer['top'],    layer['gamma_top'], layer['phi_top'], layer['Dr_top']])
            profile.append([layer['bottom'], layer['gamma_bot'], layer['phi_bot'], layer['Dr_bot']])

    z0, f_gamma, f_phi, f_Dr, f_sigma_v_eff, f_delta = sand_profile(profile)

    def qc_depth(H):
        qc0, depth0 = 0, 0
        qc1, depth1 = 8e6, 6
        qc2, depth2 = 40e6, 16
        m1 = (qc1 - qc0)/(depth1 - depth0)
        m2 = (qc2 - qc1)/(depth2 - depth1)
        if H <= depth1:
            return m1 * H
        elif H <= depth2:
            return m2 * (H - depth1) + qc1
        else:
            return qc2

    def qc_average(H, Dh):
        depths = np.linspace(max(0, H - 1.5*Dh), min(H + 1.5*Dh, 20), int(1.5*Dh/dH))
        return np.mean([qc_depth(z) for z in depths])

    def calculate_torque(Dc, Dh, Fr, ph, th, dH, H):
        gamma = f_gamma(H)
        delta_crit_rad = np.radians(delta_crit)
        a = Fr/np.tan(delta_crit_rad)
        theta = np.arctan(ph/(np.pi*Dh))
        K0 = 1 - np.sin(np.radians(32))
        Tc = np.sum([a*qc_average(z, Dh)*np.tan(delta_crit_rad)*dH*(Dc**2/2) for z in np.arange(0, H, dH)])
        Tb = qc_average(H, Dh)*np.pi*(Dc**3)*np.tan(delta_crit_rad)/12
        Th = (a*qc_average(H, Dh)*np.tan(delta_crit_rad + theta)*np.pi*(Dh**3 - Dc**3)/(12*K0) +
              a*qc_average(H, Dh)*th*np.tan(delta_crit_rad)*np.pi*(Dh**2)/2 +
              a*qc_average(H, Dh)*th*(Dh**2 - Dc**2)/8)
        return Tc + Tb + Th

    def calculate_force(Dc, Dh, Fr, ph, th, dH, H):
        gamma = f_gamma(H)
        delta_crit_rad = np.radians(delta_crit)
        a = Fr/np.tan(delta_crit_rad)
        K0 = 1 - np.sin(np.radians(32))
        Fc = np.sum([0.6*a*qc_average(z, Dh)*np.tan(delta_crit_rad)*dH*np.pi*Dc for z in np.arange(0, H, dH)])
        Fb = 0.6*qc_average(H, Dh)*np.pi*(Dc**2)/4
        Fh = (a*qc_average(H, Dh)*np.pi*(Dh**2 - Dc**2)/(4*K0) +
              a*qc_average(H, Dh)*th*np.pi*Dh/K0 +
              qc_average(H, Dh)*th*(Dh - Dc)/2)
        return Fc + Fb + Fh

    def calculate_core(T, Fy_c, Dc, tc, E, H, K=2):
        tau = 16*(T/np.pi)*(Dc/(Dc**4 - (Dc - 2*tc)**4))
        sigma_y = 4*(Fy_c/np.pi)/((Dc**2 - (Dc - 2*tc)**2))
        sigma_eq_c = np.sqrt(sigma_y**2 + 3 * tau**2)
        I = (np.pi/64)*(Dc**4 - (Dc - 2*tc)**4)
        Fy_cr = np.pi**2*E*I/(K*H)**2
        return sigma_eq_c, Fy_cr

    def calculate_helix(Fy_max, Dh, Dc, th, k):
        q = 4*Fy_max/(np.pi*(Dh**2 - Dc**2))
        return k*q*Dh**2/(4*th**2)
    
    def calculate_weld_stress(Fy_max, th, Dh, weld_length):
        M = Fy_max * th
        Q = Fy_max
        Aw = weld_length * th
        sigma_w = M / (weld_length * th**2 / 6)
        tau_w = Q / Aw
        sigma_eq_w = np.sqrt(sigma_w**2 + 3 * tau_w**2)
        return sigma_eq_w

    def calculate_axial_capacity(Dh, H):
        gamma = f_gamma(H)
        # phi_p = 6.6 + 11*np.log10(qc_average(H, Dh)/np.sqrt(gamma*H))
        phi_p = f_phi(H)
        Fps = np.tan(np.radians(psi_p)) + np.cos(np.radians(phi_p - psi_p))*\
            (np.tan(np.radians(phi_p) - np.tan(np.radians(psi_p))))
        Fs1 = 2*Fps
        Fs2 = (4/3)*Fps*np.tan(np.radians(psi_p))
        Fu = (1 + Fs1*(H/Dh) + Fs2*(H/Dh)**2)*gamma*(np.pi/4)*(Dh**2)*H
        return Fu

    Tmax = 8e6
    H_max = L
    H_values = np.linspace(0.1, H_max, 100)
    results = {
        'H': [], 'Fu': [], 'Torque': [], 'Force': [],
        'SigmaHelix': [], 'SigmaCore': [], 'BucklingLimit': [], 'SigmaWeld': [], 'FailureMode': []
    }

    for H in H_values:
        T = calculate_torque(Dc, Dh, Fr, ph, th, dH, H)
        F_inst = calculate_force(Dc, Dh, Fr, ph, th, dH, H)
        sigma_helix = calculate_helix(F_inst, Dh, Dc, th, k)
        sigma_core, Fy_cr = calculate_core(T, F_inst, Dc, tc, E, H)
        sigma_eq_weld = calculate_weld_stress(F_inst, th, Dh, weld_length=np.pi*Dc)
        Fu = calculate_axial_capacity(Dh, H)

        if T > Tmax:
            failure_mode = 'Torque limit'
        elif F_inst > Fy_cr:
            failure_mode = 'Core buckling'
        elif sigma_helix > f_y:
            failure_mode = 'Helix stress'
        elif sigma_eq_weld > f_y_weld:
            failure_mode = 'Weld stress'
        else:
            failure_mode = 'OK'

        if failure_mode == 'OK':
            results['H'].append(H)
            results['Fu'].append(Fu)
            results['Torque'].append(T)
            results['Force'].append(F_inst)
            results['SigmaHelix'].append(sigma_helix)
            results['SigmaCore'].append(sigma_core)
            results['BucklingLimit'].append(Fy_cr)
            results['SigmaWeld'].append(sigma_eq_weld)
            results['FailureMode'].append(failure_mode)

    if plot:
        plt.figure(figsize=(7, 5))
        plt.plot(results['H'], results['Fu'], label=f'Dh/Dc = {Dh/Dc:.2f}')
        plt.title("Fu vs H")
        plt.xlabel("Depth H (m)")
        plt.ylabel("Axial Capacity Fu (N)")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

    return results

if __name__ == '__main__':
    profile_map = [
        {
            'name': 'CPT_1',
            'x': 0, 'y': 0,
            'layers': [
                {'top': 0.0, 'bottom':  5.0, 'soil_type': 'sand', 'gamma_top': 10, 'gamma_bot': 11, 'phi_top': 32, 'phi_bot': 36, 'Dr_top': 55, 'Dr_bot': 65},
                {'top': 5.0, 'bottom': 20.0, 'soil_type': 'sand', 'gamma_top': 11, 'gamma_bot': 12, 'phi_top': 36, 'phi_bot': 38, 'Dr_top': 65, 'Dr_bot': 80}
            ]
        }
    ]

    results = getInstallationHelical(profile_map, 'CPT_1', D=1.5, L=10.0, d=0.5, plot=True)
      