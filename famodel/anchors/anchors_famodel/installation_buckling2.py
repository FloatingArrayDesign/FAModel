import numpy as np
import matplotlib.pyplot as plt
from support_soils import clay_profile, sand_profile
from support_pycurves import py_Matlock, py_API
from installation_suction import getInstallationSuction

def compute_Zs(s, r, t, nu):
    return (s**2/(r*t))*np.sqrt(1 - nu**2)

def compute_C(psi, rho, xi):
    return psi*np.sqrt(1 + (rho*xi/psi)**2)

def gamma_M(lam_bar):
    if lam_bar < 0.5:
        return 1.15
    elif lam_bar <= 1.0:
        return 0.85 + 0.6*lam_bar
    else:
        return 1.45

def getBucklingSuction(profile_map, location_name, D, L):
    E = 2.1e11
    fy = 325e6
    nu = 0.3
    t = (2.35 + D*20)/1e3; print(t)            # Suction pile wall thickness (m), API RP2A-WSD

    R = D/2 - t/2

    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']
    main_type = layers[0]['soil_type'].lower()
    if main_type == 'clay':
        clay_input = [[layers[0]['top'],    layers[0]['gamma_top'], layers[0]['Su_top']],
                      [layers[0]['bottom'], layers[0]['gamma_bot'], layers[0]['Su_bot']]]
        _, f_gamma, f_Su, f_sigma_v_eff, _ = clay_profile(clay_input)
    elif main_type == 'sand':
        sand_input = [[layers[0]['top'],    layers[0]['gamma_top'], layers[0]['phi_top']],
                      [layers[0]['bottom'], layers[0]['gamma_bot'], layers[0]['phi_bot']]]
        _, f_gamma, f_phi, f_Dr, f_sigma_v_eff, _ = sand_profile(sand_input)
    else:
        raise ValueError("Unsupported soil type")

    depths = np.arange(0.1, L + 0.1, 0.1)
    suction_results = getInstallationSuction(profile_map, location_name, D, L, gamma_m_install=1.5, gamma_m_retrieval=1.25)
    pe_values = np.interp(depths, suction_results['depths'], suction_results['delta_u_suction'])

    def soil_type_map(z):
        for layer in layers:
            if layer['top'] <= z <= layer['bottom']:
                return layer['soil_type'].lower()
        raise ValueError(f"No soil type defined at depth {z}")

    alpha_list = []
    y_disp = 0.001
    for z in depths:
        stype = soil_type_map(z)
        if stype == 'clay':
            py_func, _ = py_Matlock(z, D, f_gamma(z), f_Su(z), f_sigma_v_eff(z), return_curve=True)
        elif stype == 'sand':
            py_func, _ = py_API(z, D, f_phi(z), f_sigma_v_eff(z), f_Dr(z), return_curve=True)
        else:
            raise ValueError(f"Unsupported soil type at depth {z}")
        stiffness = py_func(y_disp)/y_disp if y_disp != 0 else 0
        alpha_list.append(stiffness)
    alpha_array = np.array(alpha_list)
    integral_total = np.trapz(alpha_array, depths)
    alpha_z = np.cumsum(alpha_array)*(depths[1] - depths[0])/integral_total

    UC_list = []
    PE_list = []

    for z, pe, alpha in zip(depths[:-1], pe_values[:-1], alpha_z[:-1]):
        s = L  # constant buckling length
        PE_list.append(pe)

        Zs = compute_Zs(s, R, t, nu)

        psi_a = 4.0
        xi_a = 0.702*Zs
        rho_a = 0.5*(1 + R/(150*t))**-0.5
        C_a = compute_C(psi_a, rho_a, xi_a)
        fEa = C_a*(np.pi**2*E/(12*(1 - nu**2)))*(t/s)**2

        C_m = C_a
        fEm = fEa

        psi_h = (1 + (s/L)**2)**2
        xi_h = 1.04*(s/L)*np.sqrt(Zs)
        rho_h = 0.6
        C_h = compute_C(psi_h, rho_h, xi_h)
        fEh = C_h*(np.pi**2 *E/(12*(1 - nu**2)))*(t/s)**2

        psi_t = 5.34 + (s/L)**2
        xi_t = 0.856*(s/L)*Zs**(3/4)
        rho_t = 0.6
        C_t = compute_C(psi_t, rho_t, xi_t)
        fEtau = C_t*(np.pi**2*E/(12*(1 - nu**2)))*(t/s)**2

        sigma_a = 0.5*pe*R*(1 - alpha)/t
        sigma_m = 0
        sigma_h = pe*R*(1 - alpha)/t
        tau = 0

        sigma_j = np.sqrt((sigma_a)**2 - (sigma_a)*sigma_h + sigma_h**2 + 3*tau**2)

        sigma_a0 = max(0, -sigma_a)
        sigma_m0 = max(0, -sigma_m)
        sigma_h0 = max(0, -sigma_h)

        lam_bar_sq = (fy/sigma_j)*(
            (sigma_a0/fEa) +
            (sigma_m0/fEm) +
            (sigma_h0/fEh) +
            (tau/fEtau)
        )
        lam_bar = np.sqrt(lam_bar_sq)

        gammaM = gamma_M(lam_bar)
        fksd = fy/np.sqrt(1 + lam_bar**4)/gammaM

        UC = sigma_j/fksd
        UC_list.append(UC)

    plt.figure(figsize=(6, 6))
    plt.plot(UC_list, depths[:-1], label='UC (with confinement)', color='darkred')
    plt.gca().invert_yaxis()
    plt.xlabel('Unity Check')
    plt.ylabel('Depth (m)')
    plt.title('Shell Buckling')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return {
        'depths': depths.tolist(),
        'UC_list': UC_list,
        'PE_list': PE_list,
        'alpha_z': alpha_z.tolist()
    }

if __name__ == '__main__':
    profile_map = [
        {
            'name': 'CLAY_INSTALL',
            'x': 0, 'y': 0,
            'layers': [
                {
                    'top': 0.0, 'bottom': 20.0,
                    'soil_type': 'clay',
                    'gamma_top': 9.0, 'gamma_bot': 9.0,
                    'Su_top': 5.0, 'Su_bot': 85.0
                }
            ]
        }
    ]
    getBucklingSuction(profile_map, 'CLAY_INSTALL', D=4.0, L=17.0)
