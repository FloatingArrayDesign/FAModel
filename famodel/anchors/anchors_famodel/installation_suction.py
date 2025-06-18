import numpy as np
import matplotlib.pyplot as plt
from support_soils import clay_profile

def PileWeight(Len, Dia, tw, rho):
    return ((np.pi/4)*((Dia**2 - (Dia - 2*tw)**2)*Len + (np.pi/4)*Dia**2*tw))*rho

def getInstallationSuction(profile_map, location_name, D, L, gamma_m_install=1.5, gamma_m_retrieval=1.25):
    '''
    Installation and retrieval pressure assessment for suction piles in clay using a layered profile.
    Returns a dictionary with pressure values and resistances.
    '''
    # Constants and geometry
    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']

    rhows = 66.90e3                  # Submerged steel specific weight (N/m3)
    rhow = 10e3                      # Water specific weight (N/m3) 

    WT = D/200; print(WT)
    t = (6.35 + D*20)/1e3; print(t)            # Suction pile wall thickness (m), API RP2A-WSD
    Di = D - 2*WT
    Asi = np.pi * Di
    Aso = np.pi * D
    Awall = 0.25  # m²
    Aplug = np.pi * Di**2 / 4
    Nc_strip_deep = 7.5
    Nc_circle = 9.0
    alphaD_su = 0.25
    Wp = PileWeight(L, D, WT, rhows); print(Wp)
    Wsub_steel = 750e3  # in N

    # Convert layer data into clay profile format
    z0 = layers[0]['top']
    z_bot = layers[0]['bottom']
    gamma_top = layers[0]['gamma_top']
    gamma_bot = layers[0]['gamma_bot']
    Su_top = layers[0]['Su_top']
    Su_bot = layers[0]['Su_bot']

    clay_input = [
        [z0,    gamma_top, Su_top],
        [z_bot, gamma_bot, Su_bot]
    ]

    _, f_gamma, f_Su, f_sigma_v_eff, f_alpha = clay_profile(clay_input)

    # Diagnostic: evaluate and plot alpha(z)
    z_plot = np.linspace(z0, z_bot, 100)
    alpha_plot = f_alpha(z_plot)

    # print('\n--- Adhesion Factor α(z) ---')
    # for z_check in [0, 2, 4, 6, 8, 10, 12, 14, 17]:
        # print(f"z = {z_check:>5.2f} m → α = {f_alpha(z_check):.3f}")

    # plt.figure(figsize=(5, 4))
    # plt.plot(alpha_plot, z_plot, label='α(z)', color='purple')
    # plt.gca().invert_yaxis()
    # plt.grid(True)
    # plt.xlabel('Adhesion Factor α')
    # plt.ylabel('Depth (m)')
    # plt.title('API Clay Adhesion Factor vs Depth')
    # plt.tight_layout()
    # plt.show()

    # Prepare output arrays
    depths = np.arange(0, L + 0.1, 0.1)
    Rsuction_list = []
    delta_u_suction_list = []
    delta_u_retrieval_list = []
    delta_u_all_install_list = []
    delta_u_all_retrieval_list = []
    SWP_depth = None

    for L in depths:
        z_tip = L
        z_mid = L/2
        z_tip_ext = L + alphaD_su*Di

        su_av_L = f_Su(z_mid)
        int_su = (su_av_L)*L
        su_tip = f_Su(z_tip)
        su_av_tip = f_Su(z_tip_ext)
        # alpha_i = alpha_o = float(f_alpha(z_mid))
        alpha_i = alpha_o = 0.3

        Fi = Asi*alpha_i*int_su
        Fo = Aso*alpha_o*int_su
        Qw = Awall*Nc_strip_deep*su_tip

        Rsuction = Fi + Fo + Qw
        Rretrieval = Rsuction
        delta_u_suction = max((Rsuction - Wp)/Aplug, 0.0) 
        delta_u_retrieval = (Rretrieval + Wp)/Aplug
        delta_u_all_install   = Fi/Aplug + Nc_circle*su_av_tip/gamma_m_install
        delta_u_all_retrieval = Fi/Aplug + Nc_circle*su_av_tip/gamma_m_retrieval

        Rsuction_list.append(Rsuction)
        delta_u_suction_list.append(delta_u_suction)
        delta_u_retrieval_list.append(delta_u_retrieval)
        delta_u_all_install_list.append(delta_u_all_install)
        delta_u_all_retrieval_list.append(delta_u_all_retrieval)

        if SWP_depth is None and Rsuction >= Wp:
            SWP_depth = L

    # Plotting
    fig, axs = plt.subplots(1, 3, figsize=(10, 7))

    axs[0].plot(Rsuction_list, depths, label='Installation Resistance', color='blue')
    if SWP_depth is not None:
        axs[0].axvline(Wp, color='red', linestyle='--', label=f'SWP = {SWP_depth:.2f} m')
    axs[0].set_xlabel('Installation Resistance (N)')
    axs[0].set_ylabel('Penetration (m)')
    axs[0].set_title('Installation Resistance vs Penetration')
    axs[0].grid(True)
    axs[0].invert_yaxis()
    axs[0].legend()

    axs[1].plot(delta_u_suction_list, depths, label='Underpressure', color='green')
    axs[1].plot(delta_u_all_install_list, depths, label='Δu allowable install', color='orange')
    if SWP_depth is not None:
        axs[1].axhline(SWP_depth, color='red', linestyle='--', label=f'SWP = {SWP_depth:.2f} m')
    axs[1].set_xlabel('Underpressure (Pa)')
    axs[1].set_ylabel('Penetration (m)')
    axs[1].set_title('Underpressure vs Penetration')
    axs[1].grid(True)
    axs[1].invert_yaxis()
    axs[1].legend()
    
    axs[2].plot(delta_u_retrieval_list, depths, label='Overpressure', color='green')
    axs[2].plot(delta_u_all_retrieval_list, depths, label='Δu allowable retrieve', color='orange')
    axs[2].set_xlabel('Overpressure (Pa)')
    axs[2].set_ylabel('Penetration (m)')
    axs[2].set_title('Overpressure vs Penetration')
    axs[2].grid(True)
    axs[2].invert_yaxis()
    axs[2].legend()

    plt.tight_layout()
    plt.show()

    # Final state outputs
    L = L
    z_tip = L
    z_mid = L/2
    z_tip_ext = L + alphaD_su*Di

    su_av_L = f_Su(z_mid)
    int_su = su_av_L*L
    su_tip = f_Su(z_tip)
    su_av_tip = f_Su(z_tip_ext)
    alpha_i = alpha_o = float(f_alpha(z_mid))
    alpha_i = alpha_o = 0.3

    Fi = Asi * alpha_i * int_su
    Fo = Aso * alpha_o * int_su
    Qw = Awall * Nc_strip_deep * su_tip

    Rsuction = Fi + Fo + Qw
    Rretrieval = Rsuction

    delta_u_suction = max((Rsuction - Wp)/Aplug, 0.0)  
    delta_u_retrieval = (Rretrieval + Wp)/Aplug  

    delta_u_all_install   = Fi/Aplug + Nc_circle*su_av_tip/gamma_m_install
    delta_u_all_retrieval = Fi/Aplug + Nc_circle*su_av_tip/gamma_m_retrieval

    return {
        'layers': layers,
        'depths': depths,
        'z0': z0,
        'D': D,
        'L': L,
        'Di': Di,
        'Asi': Asi,                       # m
        'Aso': Aso,                       # m
        'Aplug': Aplug,                   # m²
        'su_av_L': su_av_L/1e3,           # kPa
        'int_su': int_su/1e3,             # kN/m
        'su_tip': su_tip/1e3,             # kPa
        'su_av_tip': su_av_tip/1e3,       # kPa
        'Fi': Fi/1e3,                     # kN
        'Fo': Fo/1e3,                     # kN
        'Qw': Qw/1e3,                     # kN
        'Rsuction': Rsuction/1e6,         # MN
        'Rretrieval': Rretrieval/1e6,     # MN
        'delta_u_suction': delta_u_suction_list,              # kPa
        'delta_u_retrieval': delta_u_retrieval/1e3,           # kPa
        'delta_u_all_install': delta_u_all_install/1e3,       # kPa
        'delta_u_all_retrieval': delta_u_all_retrieval/1e3,   # kPa
        'SWP_depth': SWP_depth
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
                    'Su_top': 5.0, 'Su_bot': 45.0
                }
            ]
        }
    ]

    results = getInstallationSuction(profile_map, 'CLAY_INSTALL', D=4.0, L=17.0)
    for k, v in results.items():
        if isinstance(v, float):
            print(f"{k:<25} = {v:.2f}")
        else:
            print(f"{k:<25} = {v}")