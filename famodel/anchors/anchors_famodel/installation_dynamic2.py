import numpy as np
import matplotlib.pyplot as plt
from support_soils import clay_profile

def PileWeight(L1, L2, D1, D2, tw, gamma):
    return ((np.pi/4)*(D1**2 - (D1 - 2*tw)**2)*(L1 + L2) + 4*L2*D2*tw)*gamma

def PileVolume(L1, L2, D1, D2, tw):
    return (np.pi/4)*D1**2*(L1 + L2) + 4*L2*D2*tw

def PileWingedSurface(length, diameter1, diameter2):
    return 8*length*(diameter1 - diameter2)

def PileShaftSurface(length, diameter):
    return np.pi*diameter*length

def getInstallationDynamic(profile_map, location_name, D1, D2, L1, L2, ballast, drop_height, plot=True):
    """
    Penetration of torpedo pile using depth-based formulation (Eq. 2.10).
    """
    # Constants
    rhows = 66.90e3     # Submerged unit weight of steel (N/m³)
    rhow = 10e3         # Unit weight of water (N/m³)
    Sti = 2.0           # Installation strain-rate index (-), affects side friction strain-rate correction
    Se = 5.0            # Strain-rate multiplier (-), empirical factor for rate effects
    Ce = 0.02           # Strain-rate coefficient (-), controls shape of strain-rate correction
    Nc = 9.0            # Bearing capacity factor for undrained clay [-], used in tip resistance (q = Nc * Su)
    CD = 2.7            # Drag coefficient (-), for a cylindrical body falling in water
    dz = 1              # Depth increment (m), used in depth-stepping integration
    g = 9.81            # Gravitational acceleration (m/s²)


    # Geometry
    D = D1
    L = L1 + L2
    t = (6.35 + D2*20)/1e3

    # Soil profile
    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']
    profile = []
    for layer in layers:
        if layer['soil_type'] != 'clay': continue
        profile.append([layer['top'],    layer['gamma_top'], layer['Su_top']])
        profile.append([layer['bottom'], layer['gamma_bot'], layer['Su_bot']])
    # profile = sorted(list({(z, su, g) for z, su, g in profile}), key=lambda x: x[0])
    z0, f_gamma, f_Su, _, f_delta = clay_profile(profile)
    # print('\n--- f_gamma and f_Su vs Depth ---')
    # z_start = z0
    # z_end = profile[-1][0]
    # z_vals = np.linspace(z_start, z_end, 10)
    # for z_val in z_vals:
    #     gamma_val = f_gamma(z_val)
    #     Su_val = f_Su(z_val)
    #     print(f'z = {z_val:.1f} m → gamma = {gamma_val:.2f} N/m³, Su = {Su_val:.2f} Pa')

    # Parameters
    Af = np.pi*(D2**2)/4 + 4*(D2 - D1)*t
    As = PileWingedSurface(L1, D1, D2) + PileShaftSurface(L1 + L2, D2)
    Vol = PileVolume(L1, L2, D1, D2, t)
    Wp = PileWeight(L1, L2, D1, D2, t, rhows + rhow) + ballast; print(f'Wp = {Wp:.2f} N')
    M = Wp/g
    Mprime = M + 2*f_gamma(2*L/3)*Vol

    CD_water = 1.2
    vt = np.sqrt((2*Mprime*g)/(rhow*CD_water*Af))
    t_impact = (vt/g)*np.arccosh(np.exp(g*drop_height/vt**2))
    v_impact = vt*np.tanh(g*t_impact/vt)

    # Loop
    v = [v_impact]; #print(v)
    z = [0.0]
    term1_values, term2_values = [0.0], [0.0]
    i = 0
    while v[-1] > 0:
        zi = z[-1]
        vi = v[-1]
        Sui = f_Su(zi); #print(f'Sui = {Sui:.2f} Pa')
        gammai = f_gamma(zi); #print(f'gammai  = {gammai:.2f} N/m3')
        rhoi = gammai/g; #print(f'rhoi  = {rhoi:.2f} kg/m3')
        deltai = f_delta(zi)
        Mprime_local = M + 2*rhoi*Vol
        term1 = (Wp - Vol*gammai) - (0.5*CD*rhoi*Af*vi**2); #print(f'term1  = {term1:.2f} N')
        term2 = Sui*(Af*Nc + As*deltai/Sti)
        Se_dot = Se/(1 + (1/np.sqrt(Ce*vi/(Sui*D)) + 0.06)); #print(f'Se_dot  = {Se_dot:.2f}')
        top = 2*dz*(term1 - term2)*Se_dot; #print(f'top  = {top:.2f}')
        bottom = vi*Mprime_local; #print(f'bottom  = {bottom:.2f}')
        vi1 = vi + top/bottom
        if vi1 < 0.01:
            print(f'Stopping due to low velocity at z = {zi:.2f} m')
            break

        v.append(vi1)
        z.append(zi + dz)
        term1_values.append(term1)
        term2_values.append(term2)
        i += 1
        #print(f'z = {zi:.2f} m | term1 (drive) = {term1:.2f} N | term2 (resist) = {term2:.2f} N | net = {term1 - term2:.2f} N')

    if plot:
        fig, ax1 = plt.subplots()
        ax1.plot(z, v, 'b', label='Velocity vs Depth')
        ax1.set_xlabel('Depth (m)')
        ax1.set_ylabel('Velocity (m/s)')
        
        ax2 = ax1.twinx()
        ax2.plot(z, term1_values, 'r', label='Drive (N)')
        ax2.plot(z, term2_values, 'g', label='Resist (N)')
        ax2.set_ylabel('Forces (N)')
        
        ax1.grid(True)
        plt.title('Depth-based Penetration of Torpedo Pile')
        plt.legend(loc='upper right')
        plt.tight_layout()
        plt.show()

    return {
        'final_depth': z[-1],
        'max_velocity': max(v),
        'v_impact': v_impact,
        'steps': i
    }

if __name__ == '__main__':
    profile_map = [
        {
            'name': 'CPT_1',
            'x': 0, 'y': 0,
            'layers': [
                {'top': 0.0, 'bottom': 60.0, 'soil_type': 'clay',
                 'gamma_top': 8.5, 'gamma_bot': 8.5,
                 'Su_top': 5, 'Su_bot': 85},
                {'top': 40.0, 'bottom': 400.0, 'soil_type': 'clay',
                 'gamma_top': 8.5, 'gamma_bot': 8.5,
                 'Su_top': 85, 'Su_bot': 805}
            ]
        }
    ]

    location = 'CPT_1'
    D1 = 3.0
    D2 = 1.5
    L1 = 10.0
    L2 = 5.0
    ballast = 10000
    drop_height = 20

    results = getInstallationDynamic(profile_map, location, D1, D2, L1, L2, ballast, drop_height)

    print("\n--- Torpedo Installation Results ---")
    for k, v in results.items():
        print(f"{k}: {v:.2f}")
