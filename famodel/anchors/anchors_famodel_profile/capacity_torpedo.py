  
import numpy as np
import matplotlib.pyplot as plt
from .capacity_soils import clay_profile
from .capacity_plots import plot_torpedo

def getCapacityTorpedo(profile, soil_type, D1, D2, L1, L2, zlug, ballast, Ha, Va, plot=True):
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
    Ha : float       
        Horizontal load at pile lug elevation (N)
    Va : float          
        Vertical load at pile lug elevation (N)
    plot : bool
        Plot the capacity envelope if True

    Returns
    -------
    Dictionary with capcity, weigth and UC.
    '''

    t = (6.35 + D2*20)/1e3           # Torpedo pile wall thickness (m), API RP2A-WSD
    rhows = 66.90e3                  # Submerged steel specific weight (N/m3)
    rhow = 10e3                      # Water specific weight (N/m3) 
    
    # Average effective width
    L = L1 + L2
    A_wing_plane_1 = (D1 - D2)*L1
    A_wing_plane_2 = (D1 - D2)*np.cos(np.deg2rad(45))/2*L1
    A_shaft = D2*L
    
    # Choose based on direction:
    plane = '1'  # or '2'
    
    if plane == '1':
        Dstar = (A_wing_plane_1 + A_shaft)/L
    elif plane == '2':
        Dstar = (A_wing_plane_2 + A_shaft)/L
     
    z0, f_Su, f_sigma_v_eff, f_gamma, f_alpha = clay_profile(profile)

    a = zlug
    c = zlug + L1 + L2
    profile_depth = profile[-1, 0]

    if c > profile_depth:
        raise ValueError(
            f'Soil profile does not cover the full pile length.\n'
            f'   → Pile tip depth: {c:.2f} m\n'
            f'   → Soil profile depth: {profile_depth:.2f} m\n'
            f'Extend the soil profile to at least the pile tip depth to run the capacity model.'
        )

    z_vals = np.linspace(a, c, 100)
    Su_vals = f_Su(z_vals)
    alpha_vals = np.array([f_alpha(z) for z in z_vals])

    ez_soil = np.trapz(z_vals*Su_vals, z_vals)/np.trapz(Su_vals, z_vals)
    Su_e = f_Su(ez_soil)
    alpha_e = f_alpha(ez_soil)
    print(f"Su_e = {Su_e:.2f} kPa, ez_soil = {ez_soil:.2f} m, alpha_e = {alpha_e:.2f}")

    def PileWeight(Len1, Len2, Dia1, Dia2, tw, rho):
        return ((np.pi/4)*(Dia1**2 - (Dia1 - 2*tw)**2)*(Len1 + Len2) + 4*Len2*Dia2*tw)*rho

    def PileSurface(Len1, Len2, Dia1, Dia2):
        return np.pi*Dia1*(Len1 + Len2) + 8*Len2*Dia2*0.9

    Np_free = 3.45
    Hmax = Np_free*L*Dstar*Su_e
    Vmax = PileSurface(L1, L2, D1, D2)*alpha_e*Su_e + PileWeight(L1, L2, D1, D2, t, rhows) + ballast
    
    # Pile weight (inc. auxiliary elements) assessed as a factor
    Wp = 1.10*PileWeight(L1, L2, D1, D2, t, (rhows + rhow)) + ballast

    # Calculate actual ez_su to L ratio
    ez_ratio = (ez_soil - zlug)/L; print(f"ez_ratio = {ez_ratio:.2f} m")

    # Assign aVH and bVH based on ez_su/L
    if np.isclose(ez_ratio, 2/3, atol=0.05):
        aVH = 0.5 + L/Dstar
        bVH = 4.5 - L/(3*Dstar)
        mode = 'deep mobilization (2/3)'
    elif 0.45 <= ez_ratio <= 0.75:
        aVH = 4.5 + L/(2*Dstar)
        bVH = 3.5 - L/(4*Dstar)
        mode = 'moderate mobilization (1/2 – 3/4)'
    print(f'Interaction exponents set to aVH = {aVH:.2f}, bVH = {bVH:.2f} [{mode}]')
    
    UC = (Ha/Hmax)**aVH + (Va/Vmax)**bVH

    deg = np.linspace(0, 90, 20)
    x = np.cos(np.deg2rad(deg))
    y = (1 - x**bVH)**(1/aVH)
    X = Hmax*x; Y = Vmax*y

    if plot:
        plt.plot(X, Y, color='blue', label='VH Envelope')
        plt.plot(H, V, 'o', color='red', label='Load Point')
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
        'Unity check': UC,
        'Weight pile': Wp
    }

    return resultsTorpedo

if __name__ == '__main__':
    
    profile_clay = np.array([
        [ 0.0, 50, 8.0],
        [20.0, 50, 8.5],
        [25.0, 50, 8.5],
        [50.0, 50, 9.0]
    ])

    D1 = 3.0             # Wing diameter (m)
    D2 = 1.5             # Shaft diamter (m)
    L1 = 11.0            # Winged section length (m) 
    L2 = 10.0             # Shaft section length (m)
    zlug = 15.0          # Padeye depth (m)
    ballast = 10000      # Ballast load (N)
    H = 6.0e6            # Horizontal load (N)
    V = 8.0e6            # Vertical load (N)

    results = getCapacityTorpedo(profile_clay, 'clay', D1, D2, L1, L2, zlug, ballast, H, V, plot=True)
    print("\n--- Torpedo Pile Capacity Results ---")
    for key, val in results.items():
        print(f"{key}: {val:.2f}")
        
    plot_torpedo(profile_clay, 'clay', D1, D2, L1, L2, zlug, title='Torpedo Pile in Clay Profile')

