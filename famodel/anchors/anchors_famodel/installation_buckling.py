import numpy as np
import matplotlib.pyplot as plt
from support_soils import clay_profile, sand_profile
from installation_suction import getInstallationSuction

def compute_Zs(s, r, t, nu):
    return (s**2/(r*t))*np.sqrt(1 - nu**2)

def compute_C(psi, rho, xi):
    return psi*np.sqrt(1 + (rho*xi/psi)**2)

def gamma_M(lam_bar):
    if lam_bar < 0.5:
        return 1.15
    elif lam_bar <= 1.0:
        return 0.85 + 0.6 * lam_bar
    else:
        return 1.45

def getBucklingSuction(profile_map, location_name, D, L, t=None, fy=345e6):
    '''
    Shell buckling capacity during suction pile installation using DNV-RP-C202 and Colliard & Wallerand effective length.
    '''
    E = 2.1e11
    nu = 0.3

    if t is None:
        t = D / 200

    r = D/2 - t/2

    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']
    main_type = layers[0]['soil_type'].lower()
    if main_type == 'clay':
        clay_input = [[layers[0]['top'],    layers[0]['gamma_top'], layers[0]['Su_top']],
                      [layers[0]['bottom'], layers[0]['gamma_bot'], layers[0]['Su_bot']]]
        _, f_gamma, _, _, _ = clay_profile(clay_input)
    elif main_type == 'sand':
        sand_input = [[layers[0]['top'],    layers[0]['gamma_top'], layers[0]['phi_top']],
                      [layers[0]['bottom'], layers[0]['gamma_bot'], layers[0]['phi_bot']]]
        _, f_gamma, _, _, _ = sand_profile(sand_input)
    else:
        raise ValueError("Unsupported soil type")

    depths = np.arange(0.1, L + 0.1, 0.1)  # start from 0.1 to avoid division by zero
    suction_results = getInstallationSuction(profile_map, location_name, D, L, gamma_m_install=1.5, gamma_m_retrieval=1.25)
    pe_values = suction_results['delta_u_suction']
    pe_values = np.full_like(pe_values, 300e3)

    UC_list = []
    LB_list = []
    PE_list = []
    for z, pe in zip(depths, pe_values):
        # Colliard & Wallerand effective length
        x = L - z  # exposed length
        L_B = L*(1 + 2*(x/L) - 0.0435*(x/L)**2)
        s = L_B
        LB_list.append(s)
        PE_list.append(pe/1e3)

        Zs = compute_Zs(s, r, t, nu)

        # fEa (Axial)
        psi_a = 4.0
        xi_a = 0.702*Zs
        rho_a = 0.5*(1 + r/(150 * t))**-0.5
        C_a = compute_C(psi_a, rho_a, xi_a)
        fEa = C_a*(np.pi**2*E/(12*(1 - nu**2)))*(t/s)**2

        # fEm (Bending)
        C_m = C_a
        fEm = fEa

        # fEtau (Shear)
        psi_t = 5.34 + (s/L)**2
        xi_t = 0.856*(s/L)*Zs**(3/4)
        rho_t = 0.6
        C_t = compute_C(psi_t, rho_t, xi_t)
        fEtau = C_t*(np.pi**2*E/(12*(1 - nu**2)))*(t/s)**2
        
        # fEh (Hoop)
        psi_h = (1 + (s/L)**2)**2
        xi_h = 1.04*(s/L)*np.sqrt(Zs)
        rho_h = 0.6
        C_h = compute_C(psi_h, rho_h, xi_h)
        fEh = C_h*(np.pi**2*E/(12*(1 - nu**2)))*(t/s)**2

        sigma_a = 0.5*pe*r/t
        sigma_m = 0
        sigma_h = pe*r/t
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

    fig, axs = plt.subplots(1, 3, figsize=(18, 6))

    axs[0].plot(UC_list, depths, label='UC', color='blue')
    axs[0].invert_yaxis()
    axs[0].set_xlabel('Unity check')
    axs[0].set_ylabel('Depth (m)')
    axs[0].set_title('Shell buckling UC vs. Depth')
    axs[0].grid(True)
    axs[0].legend()

    axs[1].plot(LB_list, depths, label='Buckling Length (L_B)', color='blue')
    axs[1].invert_yaxis()
    axs[1].set_xlabel('Effective buckling length (m)')
    axs[1].set_ylabel('Depth (m)')
    axs[1].set_title('Buckling Length vs. Depth')
    axs[1].grid(True)
    axs[1].legend()
    
    axs[2].plot(PE_list, depths, label='Underpressure', color='green')
    axs[2].invert_yaxis()
    axs[2].set_xlabel('Underpressure (kPa)')
    axs[2].set_ylabel('Depth (m)')
    axs[2].set_title('Suction Pressure vs. Depth')
    axs[2].grid(True)
    axs[2].legend()

    plt.tight_layout()
    plt.show()

    return {
        'depths': depths.tolist(),
        'UC_list': UC_list,
        'LB_list': LB_list
    }

if __name__ == '__main__':
    profile_map = [
        {
            'name': 'CPT_1',
            'x': 0, 'y': 0,
            'layers': [
                {
                    'top': 0.0, 'bottom': 20.0,
                    'soil_type': 'clay',
                    'gamma_top': 9.0, 'gamma_bot': 9.0,
                    'Su_top': 5.0, 'Su_bot': 45.0
                }
            ]
        }
    ]
    getBucklingSuction(profile_map, 'CPT_1', D=2.0, L=10.4)
