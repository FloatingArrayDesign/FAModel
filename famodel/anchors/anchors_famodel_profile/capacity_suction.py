   
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
from .capacity_soils import clay_profile, sand_profile
from .capacity_plots import plot_suction

def getCapacitySuction(profile, soil_type, D, L, zlug, Ha, Va, plot=True):    
    '''Calculate the inclined load capacity of a suction pile in sand or clay following S. Kay methodology.
    The calculation is based on the soil profile, anchor geometry and inclined load.  
    
    Parameters
    ----------
    profile : array
        Soil profile as a 2D array: (z, parameters)
            Clay soil profile (z (m), Su (kPa), gamma (kN/m³))
            Sand soil profile (z (m), phi (deg), gamma (kN/m³), Dr (%))
    soil_type : string
        Select soil condition, 'clay' or 'sand'
    D : float 
        Suction pile diameter (m)
    L : float 
        Suction pile length from pile head (m)
    zlug: float
        Embedded depth of the main padeye (m)
    thetalug: float
        Angle of tilt misaligment (deg) (default value: 5.0)
    psilug: float
        Angle of twist misaligment (deg) (default value: 7.5)
    Ha : float       
        Horizontal load at pile lug elevation (N)
    Va : float          
        Vertical load at pile lug elevation (N)
    plot : bool
        Plot the capacity envelope if True
                    
    Returns
    -------
    Dictionary with capcity, weigths and UC.
    ''' 
    
    z0 = profile[0][0]      
    lambdap = (L - z0)/D; m = 2/3;   # Suction pile slenderness ratio
    t = (6.35 + D*20)/1e3            # Suction pile wall thickness (m), API RP2A-WSD
    rlug = D/2                       # Radial position of the lug
    thetalug = 5                     # Angle of tilt misaligment, default is 5. (deg)
    psilug = 7.5                     # Angle of twist misaligment, default is 7.5. (deg)
    rhows = 66.90e3                  # Submerged steel specific weight (N/m3)
    rhow = 10e3                      # Water specific weight (N/m3) 
    
    # Outer and inner surface of the pile skirt
    def PileSurface(Len, Dia):
        Sp = np.pi*Dia*Len
        return Sp    
    # Dry and wet mass of the pile    
    def PileWeight(Len, Dia, tw, rho):
        Wp = ((np.pi/4)*((Dia**2 - (Dia - 2*tw)**2)*Len + (np.pi/4)*Dia**2*tw))*rho
        return Wp 
    # Mass of the soil plug      
    def SoilWeight(Len, Dia, tw, gamma_soil): 
        Wsoil =(np.pi/4)*(Dia - 2*tw)**2*Len*gamma_soil
        return Wsoil
    # Tilt and twist effects due to installation misaligments
    def rlugTilt(r, z, theta):
        R = r*np.cos(np.deg2rad(theta)) - z*np.sin(np.deg2rad(theta))
        return R   
    def zlugTilt(r, z, theta):
        Z = r*np.sin(np.deg2rad(theta)) + z*np.cos(np.deg2rad(theta))
        return Z    
    
    if soil_type == 'clay':
        z0, f_Su, f_sigma_v_eff, f_gamma, f_alpha = clay_profile(profile)
        z_vals = np.linspace(z0, L, 10)

        Su_vals = f_Su(z_vals)
        ez = np.trapz(z_vals*Su_vals, z_vals)/np.trapz(Su_vals, z_vals); 
        ez_soil = ez; print(f"ez_soil = {ez_soil:.2f} m")
        gamma_vals = f_gamma(z_vals)

        Su_av_L = f_Su(ez_soil)
        Su_tip = f_Su(L); # print(f"Su_tip = {Su_tip:.2f} Pa")
        sigma_v_eff = f_sigma_v_eff(ez_soil)
        gamma_av = np.trapz(gamma_vals, z_vals)/(L - z0); # print(f"gamma_av = {gamma_av:.2f} kN/m³")
        alpha_av = float(f_alpha(ez_soil))
                
        Nc = min(6.2*(1 + 0.34*np.arctan(lambdap)), 9)
        Np_fixed = 10.25; Np_free = 4 
        Hmax = Np_fixed*(L - z0)*D*Su_av_L; print(f'Su_av_L = {Su_av_L:.2f} Pa')
        print(f'Hmax = {Hmax:.2f} N')
        
        M = -Va*rlugTilt(rlug, zlug, thetalug) - Ha*(zlugTilt(rlug, zlug, thetalug) - ez_soil)
        # print(f"rlug_eff = {rlugTilt(rlug, zlug, thetalug):.2f} m")
        # print(f"zlug_eff = {zlugTilt(rlug, zlug, thetalug):.2f} m")
        print(f"M = {M:.2f} Nm")

        # --- MH Ellipse Parameters for Clay (Kay 2014) ---
        # ΔφMH (piecewise based on L/D)
        if 0.5 <= lambdap < 1.125:
            delta_phi = 0.32 + 4.32*lambdap; # print(delta_phi)
        elif 1.125 <= lambdap < 2.0:
            delta_phi = 7.13 - 1.71*lambdap; # print(delta_phi)
        elif 2.0 <= lambdap <= 8.0:
            delta_phi = 4.55 - 0.425*lambdap; # print(delta_phi)
        else:
            raise ValueError('L/D ratio out of bounds for MH ellipse formulation.')

        phi_MH = -np.arctan(ez_soil/(L -z0)) - np.deg2rad(delta_phi)
        a_MH = Np_fixed/np.cos(phi_MH)
        delta_bMH = 0.45*(lambdap)**(-0.9) if lambdap <= 1.5 else 0
        b_MH = -Np_free*np.sin(phi_MH) + delta_bMH
        print('M cos(phi)/a_MH =', (M*np.cos(phi_MH))/a_MH)
        print('M sin(phi)/b_MH =', (M*np.sin(phi_MH))/b_MH)

        # Solve MH ellipse for Hmax
        def f(Hmax):
            term1 = ((M*np.cos(phi_MH) + Hmax*np.sin(phi_MH))/a_MH)**2
            term2 = ((M*np.sin(phi_MH) - Hmax*np.cos(phi_MH))/b_MH)**2
            return term1 + term2 - 1

        try:
            Hmax = max(fsolve(f, Hmax*0.8)[0], 0.0)
        except:
            Hmax = 0.0
            
        print(f'Hmax (MH ellipse) = {Hmax:.2f} N')

        To = PileSurface((L - z0), D)*alpha_av*Su_av_L
        Ti = PileSurface((L - z0), D - 2*t)*alpha_av*Su_av_L
        Tbase = np.pi*D**3*Su_tip/12
        Tmax = min(To + Ti, To + Tbase)

        T = Ha*rlug*np.sin(np.deg2rad(psilug))
        nhuT = T/Tmax
        nhuV = Ha/To
        nhuVstar = np.sqrt(nhuV**2 - nhuT**2)
        alphastar = alpha_av*(nhuVstar/nhuV); print(f"alphastar   = {alphastar:.3f}")
        
        Vmax1 = PileWeight(L, D, t, rhows) + PileSurface((L - z0), D)*alphastar*Su_av_L + Nc*Su_tip*(np.pi/4)*D**2
        Vmax2 = PileWeight(L, D, t, rhows) + PileSurface((L - z0), D)*alphastar*Su_av_L + PileSurface((L - z0), D - 2*t)*alphastar*Su_av_L
        Vmax3 = PileWeight(L, D, t, rhows) + PileSurface((L - z0), D)*alphastar*Su_av_L + SoilWeight((L - z0), D, t, gamma_av)
        Vmax = min(Vmax1, Vmax2, Vmax3)
        
        print(f"Vmax1    = {Vmax1:.2f} N"); print(f"Vmax2    = {Vmax2:.2f} N"); print(f"Vmax3    = {Vmax3:.2f} N")
        

    elif soil_type == 'sand':
        z0, f_phi, f_sigma_v_eff, f_gamma, f_Dr, f_delta = sand_profile(profile)
        z_vals = np.linspace(z0, L, 10)

        sigma_v_vals = f_sigma_v_eff(z_vals)
        ez = np.trapz(z_vals*sigma_v_vals, z_vals) / np.trapz(sigma_v_vals, z_vals)
        ez_soil = ez - z0        

        phi_vals = f_phi(z_vals)
        gamma_vals = f_gamma(z_vals)
        sigma_v_vals = f_sigma_v_eff(z_vals)
        Dr_vals = f_Dr(z_vals)
        delta_vals = f_delta(z_vals)

        phi_av = np.trapz(phi_vals, z_vals)/(L - z0)
        gamma_av = np.trapz(gamma_vals, z_vals)/(L - z0)
        delta_av = np.trapz(delta_vals, z_vals)/(L - z0)
        sigma_av_L = np.trapz(sigma_v_vals, z_vals)/(L - z0)
        sigma_tip = f_sigma_v_eff(L)

        Nq = np.e**(np.pi*np.tan(np.radians(phi_av)))*(np.tan(np.radians(45) + np.radians(phi_av)/2))**2
        Hmax = 0.5*D*Nq*gamma_av*(L - z0)**2
        Np_free = 3.0  

        M = -Va*rlugTilt(rlug, zlug, thetalug) - Ha*(zlugTilt(rlug, zlug, thetalug) - ez_soil)

        # --- MH Ellipse Parameters for Clay (Kay 2014) ---
        # ΔφMH (piecewise based on L/D)
        if 0.5 <= lambdap < 1.125:
            delta_phi = 0.32 + 4.32*lambdap
        elif 1.125 <= lambdap < 2.0:
            delta_phi = 7.13 - 1.71*lambdap
        elif 2.0 <= lambdap <= 6.0:
            delta_phi = 4.55 - 0.425*lambdap
        else:
            raise ValueError('L/D ratio out of bounds for MH ellipse formulation.')

        phi_MH = -np.arctan(ez_soil/(L - z0)) - np.deg2rad(delta_phi)
        a_MH = Nq/np.cos(phi_MH)     
        delta_bMH = 0.45*(lambdap)**(-0.9) if lambdap <= 1.5 else 0
        b_MH = -Nq*np.sin(phi_MH) + delta_bMH

        # Solve MH ellipse for Hmax
        def f(Hmax):
            term1 = ((M*np.cos(phi_MH) + Hmax*np.sin(phi_MH))/a_MH)**2
            term2 = ((M*np.sin(phi_MH) - Hmax*np.cos(phi_MH))/b_MH)**2
            return term1 + term2 - 1

        Hmax = fsolve(f, 0.8*Hmax)[0]
        print(f'Hmax (MH ellipse) = {Hmax:.2f} N')

        To = PileSurface((L - z0), D)*delta_av*sigma_av_L
        Ti = PileSurface((L - z0), D - 2*t)*delta_av*sigma_av_L
        Tbase = np.pi*D**3*sigma_tip/12
        Tmax = min(To + Ti, To + Tbase)

        T = Ha*rlug*np.sin(np.deg2rad(psilug))
        Fo = delta_av*sigma_av_L*L*np.pi*D
        nhuT = T/Tmax
        nhuV = Ha/Fo
        nhuVstar = np.sqrt(nhuV**2 - nhuT**2)
        deltastar = delta_av*(nhuVstar/nhuV)

        Vmax2 = PileWeight(L, D, t, rhows) + PileSurface((L - z0), D)*deltastar*sigma_av_L + PileSurface((L - z0), D - 2*t)*deltastar*sigma_av_L
        Vmax3 = PileWeight(L, D, t, rhows) + PileSurface((L - z0), D)*deltastar*sigma_av_L + SoilWeight((L - z0), D, t, gamma_av)
        Vmax = min(Vmax2, Vmax3)
        
    # Pile weight (inc. stiffening plus vent) assessed as a factor
    Wp = 1.10*PileWeight(L, D, t, (rhows + rhow)) 
    # Submerged weight of the soil plug
    Wsoil = SoilWeight((L - z0), D, t, gamma_av)
    
    # Capacity envelope
    aVH = 0.5 + lambdap; bVH = 4.5 + lambdap/3 
    # print('Env. exp = ' +str(aVH)+'   '+str(bVH))
    UC = (Ha/Hmax)**aVH + (Va/Vmax)**bVH      
    x = np.cos(np.linspace (0, np.pi/2, 100))
    y = (1 - x**bVH)**(1/aVH)
    X = Hmax*x; Y = Vmax*y
    
    if plot:
        plt.figure(figsize=(6, 5))
        plt.plot(X, Y, color = 'b', label='VH Envelope')
        plt.plot(Ha, Va, 'o', color = 'r', label='Load Point')
            
        # Set labels and title
        plt.xlabel('Horizontal capacity (N)')
        plt.ylabel('Vertical capacity (N)')
        plt.suptitle('VH suction pile capacity envelope')
        plt.axis([0, 1.3*max(X[0], Ha), 0, 1.3*max(Y[-1], Va)]) 
        plt.legend()
        plt.grid(True)
        plt.show()
    
    resultsSuction = {
        'Horizontal max.': Hmax,
        'Vertical max.': Vmax,
        'Unity check': UC,
        'Weight pile': Wp,
        'Weight soil': Wsoil,
        't': t
    }
    
    return resultsSuction
    
if __name__ == '__main__':

    # Clay profile: [depth (m), Su (kPa), gamma (kN/m³)]
    profile_clay = np.array([
        [ 2.0,  25,  8.5],
        [ 4.0,  25,  8.5],
        [ 6.0,  25,  8.5],
        [20.0,  25,  8.5]
    ])
    
    # Sand profile: [depth (m), phi (deg), gamma (kN/m³), Dr(%)]
    profile_sand = np.array([
        [ 1.0, 28, 8.0, 75],
        [ 5.0, 35, 8.5, 75],
        [ 8.0, 38, 9.0, 75],
        [20.0, 42, 9.5, 75]
    ])

    D = 3.0             # Diameter (m)
    L = 10.0            # Length (m)
    zlug = 8.0          # Padeye depth (m)
    Ha = 2.0e6          # Horizontal load (N)
    Va = 3.0e6          # Vertical load (N)

    # === CLAY ===
    resultsSuction_clay = getCapacitySuction(profile_clay, 'clay', D, L, zlug, Ha, Va, plot=True)
    for key, val in resultsSuction_clay.items(): 
        print(f"{key}: {val:.2f}")
        
    # Plot suction pile with the clay profile
    profile_clay_plot = [(float(z), float(Su), 'clay') for z, Su, _ in profile_clay]
    plot_suction(profile_clay_plot, 'clay', L, D, zlug)
      
    # === SAND ===
    # resultsSuction_sand = getCapacitySuction(profile_sand, 'sand', D, L, zlug, Ha, Va, plot=True)
    # for key, val in resultsSuction_sand.items(): 
    #     print(f"{key}: {val:.2f}") 
    
    # # Sand profile formatted for plotting
    # profile_sand_plot = [(float(z), float(phi), 'sand') for z, phi, _, _ in profile_sand]
    # plot_suction(profile_sand_plot, 'sand', L, D, zlug, title='Suction Pile in Sand Profile')

