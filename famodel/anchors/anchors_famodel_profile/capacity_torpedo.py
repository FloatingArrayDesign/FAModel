  
import numpy as np
import matplotlib.pyplot as plt
from capacity_soils import clay_profile
from capacity_plots import plot_torpedo

def getCapacityTorpedo(profile, soil_type, D1, D2, L1, L2, zlug, ballast, H, V, plot=True):
    '''Calculate the inclined load capacity of a torpedo pile in clay following S. Kay methodology.
    The calculation is based on the soil profile, anchor geometry and inclined load.  

    Parameters
    ----------
    profile : array
        Clay soil profile (z, Su, gamma)
            Clay soil profile (z (m), Su (kPa), gamma (kN/m³))
    soil_type : string
        Select soil condition, 'clay' 
    D1 : float 
        Wing diameter (m)
    D2 : float    
        Shaft diameter (m)
    L1 : float 
        Winged section length (m)
    L2 : float 
        Shaft section length (m)
    zlug : float 
        Padeye embedment depth (m)       
    H : float
        Horizontal load at padeye depth (N)
    V : float
        Vertical load at padeye depth (N)
    plot : bool
        Plot the capacity envelope if True

    Returns
    -------
    Dictionary with capcity, weigth and UC.
    '''

    t = (6.35 + D2*20)/1e3           # Torpedo pile wall thickness (m), API RP2A-WSD
    rhows = 66.90e3                  # Submerged steel specific weight (N/m3)
    rhow = 10e3                      # Water specific weight (N/m3) 
    
    L = L1 + L2
    Dstar = (D1*L1 + (D1 + 2*D2)*L2)/L
    lambdap = L/Dstar

    z0, f_Su, f_sigma_v_eff, f_gamma, f_alpha = clay_profile(profile)

    a = zlug
    b = zlug + L1
    c = zlug + L1 + L2

    z_vals = np.linspace(a, min(c, profile[-1, 0]), 200)
    Su_vals = f_Su(z_vals)
    alpha_vals = np.array([f_alpha(z) for z in z_vals])

    ez_Su = np.trapz(z_vals * Su_vals, z_vals)/np.trapz(Su_vals, z_vals)
    z_target = min(zlug + ez_Su, profile[-1, 0])
    Su_e = f_Su(z_target)
    alpha_e = f_alpha(z_target)

    def PileWeight(Len1, Len2, Dia1, Dia2, tw, rho):
        return ((np.pi/4)*(Dia1**2 - (Dia1 - 2*tw)**2)*(Len1 + Len2) + 4*Len2*Dia2*tw)*rho

    def PileSurface(Len1, Len2, Dia1, Dia2):
        return np.pi*Dia1*(Len1 + Len2) + 8*Len2*Dia2*0.9

    Np_free = 3.45
    Hmax = L*Dstar*Np_free*Su_e
    Vmax = PileSurface(L1, L2, D1, D2)*alpha_e*Su_e + PileWeight(L1, L2, D1, D2, t, rhows) + ballast
    
    # Pile weight (inc. auxiliary elements) assessed as a factor
    Wp = 1.10*PileWeight(L1, L2, D1, D2, t, (rhows + rhow)) 

    aVH = 4.5 + L/(2*Dstar)
    bVH = 3.5 - L/(4*Dstar)
    UC = (H/Hmax)**aVH + (V/Vmax)**bVH

    deg = np.linspace(0, 90, 20)
    x = np.cos(np.deg2rad(deg))
    y = (1 - x**bVH)**(1/aVH)
    X = Hmax*x; Y = Vmax*y

    if plot:
        plt.plot(X, Y, color='blue', label='VH Envelope')
        plt.plot(H, V, '*', color='red', label='Load Point')
        plt.xlabel('Horizontal Load (N)')
        plt.ylabel('Vertical Load (N)')
        plt.title('VH torpedo pile capacity envelope')
        plt.grid(True)
        plt.legend()
        plt.axis([0, 1.3*max(X[0], H), 0, 1.3*max(Y[-1], V)])
        plt.show()

    resultsTorpedo = {
        'Horizontal max.': Hmax,
        'Vertical max.': Vmax,
        'Unity Check': UC,
        'Weight pile': Wp
    }

    return resultsTorpedo

if __name__ == '__main__':
    
    profile_clay = np.array([
        [ 0.0, 10, 8.0],
        [20.0, 25, 8.5],
        [28.0, 45, 8.5],
        [50.0, 50, 9.0]
    ])

    D1 = 2.5             # Wing diameter (m)
    D2 = 0.8             # Shaft diamter (m)
    L1 = 15.0            # Winged section length (m) 
    L2 = 5.0             # Shaft section length (m)
    zlug = 18.0          # Padeye depth (m)
    ballast = 10000      # Ballast load (N)
    H = 6.0e6            # Horizontal load (N)
    V = 8.0e6            # Vertical load (N)

    results = getCapacityTorpedo(profile_clay, 'clay', D1, D2, L1, L2, zlug, ballast, H, V, plot=True)
    print("\n--- Torpedo Pile Capacity Results ---")
    for key, val in results.items():
        print(f"{key}: {val:.2f}")
        
    plot_torpedo(profile_clay, 'clay', D1, D2, L1, L2, zlug, title='Torpedo Pile in Clay Profile')

