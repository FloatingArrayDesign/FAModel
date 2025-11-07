
import numpy as np
import matplotlib.pyplot as plt
from support_soils import clay_profile

def PileWeight(L1, L2, D1, D2, tw, rho):
    return ((np.pi/4)*(D1**2 - (D1 - 2*tw)**2)*(L1 + L2) + 4*L2*D2*tw)*rho

def PileVolume(L1, L2, D1, D2, tw):
    return (np.pi/4)*D1**2*(L1 + L2) + 4*L2*D2*tw

def PileWingedSurface(length, diameter1, diameter2):
    return 8*length*(diameter1 - diameter2)

def PileShaftSurface(length, diameter):
    return np.pi*diameter*length

def getInstallationDynamic(profile_map, location_name, D1, D2, L1, L2, ballast, drop_height, plot=True):
    """
    Deterministic installation model of a torpedo pile in clay based on time-domain integration.
    Implements the model by True (1976) as adapted in Kazue et al. (2020), accounting for layered soil.
    """
    # Constants
    rhows = 66.90e3                  # Submerged steel specific weight (N/m3)
    rhow = 10e3                      # Water specific weight (N/m3) 
    Sti = 3.0
    Se = 5.0
    Ce = 0.02
    delta = 0.9
    Nc = 9.0
    CD = 2.7
    dt = 0.002
    tmax = 15.0
    beta = 28/27
    g = 9.81

    # Geometry
    D = D1                  # use the wing diameter for frontal area
    L = L1 + L2
    t = (6.35 + D2*20)/1e3  # assumed same as in getCapacityTorpedo

    # Retrieve soil profile and construct full profile list
    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']
    profile = []
    for layer in layers:
        if layer['soil_type'] != 'clay':
            continue
        profile.append([layer['top'],    layer['Su_top'], layer['gamma_top']])
        profile.append([layer['bottom'], layer['Su_bot'], layer['gamma_bot']])

    # Sort and remove duplicates if needed
    # profile = sorted(list({(z, su, g) for z, su, g in profile}), key=lambda x: x[0])
    z_ref, f_gamma, f_Su, _, f_delta = clay_profile(profile)

    # Precompute parameters
    Af = np.pi*(D2**2)/4 + 4*(D2 - D1)*t
    As = PileWingedSurface(L1, D1, D2) + PileShaftSurface(L1 + L2, D2)
    Vol = PileVolume(L1, L2, D1, D2, t)
    Wp = PileWeight(L1, L2, D1, D2, t, rhows) + ballast; #print(f'Wp = {Wp:.2f} N')
    M = Wp/g
    Mprime = M + 2*rhow*Vol

    # Closed-form solution for v_impact from free fall in water
    CD_water = 1.2
    A_water = Af
    vt = np.sqrt((2*Mprime*g)/(rhow*CD_water*A_water))
    t_impact = (vt/g)*np.arccosh(np.exp(g*drop_height/vt**2))
    v_impact = vt*np.tanh(g*t_impact/vt)

    # Initial conditions
    t = [0.0]
    z = [0.0]
    v = [v_impact]
    a = [Wp/Mprime]

    # Integration constants
    beta1 = dt**2*(0.5 - beta)
    beta2 = dt**2*beta
    gamma1 = -0.5*dt
    gamma2 = 1.5*dt
    nsteps = int(tmax/dt)

    for i in range(nsteps):
        zn = z[-1]
        vn = v[-1]
        an = a[-1]

        Su = f_Su(zn)
        gamma = f_gamma(zn)
        delta = f_delta(zn)
        rho = gamma/g
        Se_dot = Se/(1 + (1/np.sqrt(Ce*vn/(Su*D2)) + 0.06)) if vn > 0 else 1.0

        Mprime_local = M + 2*rho*Vol
        FD = 0.5*rhow*CD*Af*vn*abs(vn)
        FT = Su*Nc*Af*Se_dot
        FS = Su*As*delta*Se_dot/Sti

        f_total = (Wp - Vol*gamma) - FD - FT - FS 
        an1 = f_total/Mprime_local

        zn1 = zn + dt*vn + beta1*an + beta2*an1
        vn1 = vn + gamma1*an + gamma2*an1

        if vn1 < 0:
            break  # penetration stops

        t.append(t[-1] + dt)
        z.append(zn1)
        v.append(vn1)
        a.append(an1)

    if plot:
        fig, ax1 = plt.subplots()
        ax1.plot(t, z, 'b', label='Penetration depth (m)')
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Depth (m)')
        ax1.grid(True)

        ax2 = ax1.twinx()
        ax2.plot(t, v, 'r', label='Velocity (m/s)')
        ax2.set_ylabel('Velocity (m/s)')

        fig.suptitle('Torpedo Pile Installation Response')
        fig.legend(loc='upper right')
        plt.tight_layout()
        plt.show()

    return {
        'final_depth': z[-1],
        'max_velocity': max(v),
        'penetration_time': t[-1]
    }

if __name__ == '__main__':
    profile_map = [
        {
            'name': 'CPT_1',
            'x': 0, 'y': 0,
            'layers': [
                {'top': 0.0, 'bottom': 30.0, 'soil_type': 'clay',
                 'gamma_top': 8.0, 'gamma_bot': 8.5,
                 'Su_top': 5, 'Su_bot': 20},
                {'top': 30.0, 'bottom': 100.0, 'soil_type': 'clay',
                 'gamma_top': 8.5, 'gamma_bot': 9.0,
                 'Su_top': 20, 'Su_bot': 60}
            ]
        }
    ]

    location = 'CPT_1'
    D1 = 3.0
    D2 = 1.5
    L1 = 5.0
    L2 = 5.0
    ballast = 350000
    drop_height = 200

    results = getInstallationDynamic(profile_map, location, D1, D2, L1, L2, ballast, drop_height)

    print("\n--- Torpedo Installation Results ---")
    for k, v in results.items():
        print(f"{k}: {v:.2f}")
