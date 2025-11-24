
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from .support_soils import clay_profile, sand_profile

def getCapacitySuction(profile_map, location_name, D, L, zlug, Ha, Va, thetalug=5, psilug=7.5, plot=False, display=0):
    '''Calculate the inclined load capacity of a suction pile in sand or clay following S. Kay methodology.
    The calculation is based on the soil profile, anchor geometry and inclined load.  
    
    Parameters
    ----------
    profile : array
        Soil profile as a 2D array: (z, parameters)
            Clay soil profile (z (m), Su (kPa), gamma (kN/m³))
            Sand soil profile (z (m), phi (deg), gamma (kN/m³), Dr (%))
    location_name : str
        Name of the location in profile_map (e.g. 'CPT_1')
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
    
    # Retrieve soil layers from map
    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']
    
    z0 = layers[0]['top']            # Mudline elevation    
    lambdap = (L - z0)/D             # Suction pile slenderness ratio
    t = (6.35 + D*20)/1e3            # Suction pile wall thickness (m), API RP2A-WSD
    rlug = D/2                       # Radial position of the lug
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
    # Ellipse crossing with constant values
    def horizontal_cross(H, M, M_target):
        crossings = []
        for i in range(len(M_rot) - 1):
            if (M[i] - M_target) * (M[i+1] - M_target) < 0:
                # Interpolation to get more precise value at crossing
                H_cross = np.interp(M_target, [M[i], M[i+1]], [H[i], H[i+1]])
                crossings.append(H_cross)
        return crossings  
    def vertical_cross(H, M, H_target):
        crossings = []
        for i in range(len(H) - 1):
            if (H[i] - H_target) * (H[i+1] - H_target) < 0:
                # Interpolation to get more precise value at crossing
                M_cross = np.interp(H_target, [H[i], H[i+1]], [M[i], M[i+1]])
                crossings.append(M_cross)
        return crossings
    
    Np_fixed = 11.65
    Np_free = 3.5
    Nc = min(6.2*(1 + 0.34*np.arctan(lambdap)), 9)
    
    # Initialize
    sum_ez_weighted = 0.0
    Hmax_final = 0.0
    Vmax_final = 0.0
    layer_data = []

    # Profile check points
    npts = 10

    for layer in layers:
        soil_type = layer['soil_type']
        z_top = layer['top']
        z_bot = layer['bottom']

        if soil_type == 'clay':
            # Prepare soil profile for clay
            profile = [
                [z_top, layer['gamma_top'], layer['Su_top']],
                [z_bot, layer['gamma_bot'], layer['Su_bot']]
            ]
    
            z_ref, f_gamma, f_Su, f_sigma_v_eff, f_alpha = clay_profile(profile)
    
            # Clip the layer first
            z_top_clip = max(z_top, z0)
            z_bot_clip = min(z_bot, z0 + (L - z0))
            dz_clip = z_bot_clip - z_top_clip
            if display > 1: print(f'dz_clip  = {dz_clip:.2f} m')
    
            if dz_clip <= 0:
                continue  # Skip layers fully above or below
    
            # Calculate properties over clipped dz
            z_vals = np.linspace(z_top_clip, z_bot_clip, npts)
            Su_vals = f_Su(z_vals)
            Su_total = np.trapz(Su_vals, z_vals)
            Su_moment = np.trapz(Su_vals*z_vals, z_vals)
            
            ez_layer = Su_moment/Su_total
            Su_av_z = f_Su(ez_layer)
            
            if display > 1: print(f'ez_layer = {ez_layer:.2f} m')
            if display > 1: print(f'Su_av_z (at ez_layer) = {Su_av_z:.2f} Pa')
            
            Su_bot = f_Su(z_bot_clip)
            gamma_vals = f_gamma(z_vals)
            gamma_av = np.mean(gamma_vals)

            # Calculate Hmax for clay
            Hmax_layer = Np_fixed*D*dz_clip*Su_av_z
    
            layer_data.append({
                'z_top': z_top_clip,
                'z_bot': z_bot_clip,
                'dz': dz_clip,
                'Hmax_layer': Hmax_layer,
                'ez_layer': ez_layer
            })
    
            sigma_v_eff = f_sigma_v_eff(np.mean(z_vals))
            alpha_av = float(f_alpha(np.mean(z_vals)))
            
            # Side shear To and Ti
            To = PileSurface(dz_clip, D)*alpha_av*Su_av_z
            Ti = PileSurface(dz_clip, D - 2*t)*alpha_av*Su_av_z
            
            # Tip resistance
            if abs(z_bot_clip - (z0 + (L - z0))) < 1e-3:  # tip check
                Tbase = (np.pi/12)*D**3*Su_bot
            else:
                Tbase = 0.0
            
            Tmax = min(To + Ti, To + Tbase)
            
            # Torque induced by horizontal load
            T = Ha*rlug*np.sin(np.deg2rad(psilug))
            
            nhuT = T/Tmax 
            nhuV = Ha/To 
            nhuVstar = np.sqrt(nhuV**2 - nhuT**2) 
            alphastar = alpha_av*(nhuVstar/nhuV)
            if display > 1: print(f"alphastar   = {alphastar:.3f}")
            
            # Constant weight
            Pile_Head = PileWeight(z0, D, t, rhows)
            
            # Vertical failure modes
            Vmax1 = None
            if np.isclose(z_bot_clip, z0 + (L - z0), atol=0.1):
                Vmax1 = PileWeight(dz_clip, D, t, rhows) + PileSurface(dz_clip, D)*alphastar*Su_av_z + Nc*Su_bot*(np.pi/4)*D**2
            # else:
            #     Vmax1 = np.inf  # No tip resistance unless at tip
            
            Vmax2 = PileWeight(dz_clip, D, t, rhows) + PileSurface(dz_clip, D)*alphastar*Su_av_z + PileSurface(dz_clip, D - 2*t)*alphastar*Su_av_z
            Vmax3 = PileWeight(dz_clip, D, t, rhows) + PileSurface(dz_clip, D)*alphastar*Su_av_z + SoilWeight(dz_clip, D, t, gamma_av)
            
            Vmax_candidates = [v for v in [Vmax1, Vmax2, Vmax3] if v is not None]
            Vmax_layer = max(Vmax_candidates)

            # Sum vertical capacities
            Vmax_final += Vmax_layer
            
            # Print layer debug info
            if display > 0: print(f"Vmax_layer    = {Vmax_layer:.2f} N")
            if display > 0: print(f"Vmax1   = {Vmax1:.2f} N" if Vmax1 is not None else "Vmax1   = not applicable")
            if display > 0: print(f"Vmax2   = {Vmax2:.2f} N")
            if display > 0: print(f"Vmax3   = {Vmax3:.2f} N")
            
        elif soil_type == 'sand':
            # Prepare soil profile for sand
            profile = [
                [z_top, layer['gamma_top'], layer['phi_top'], layer['Dr_top']],
                [z_bot, layer['gamma_bot'], layer['phi_bot'], layer['Dr_bot']]
            ]
        
            z_ref, f_gamma, f_phi, _, f_sigma_v_eff, f_delta = sand_profile(profile)
        
            # Clip the layer within pile embedded length
            z_top_clip = max(z_top, z0)
            z_bot_clip = min(z_bot, z0 + (L - z0))
            dz_clip = z_bot_clip - z_top_clip
        
            if dz_clip <= 0:
                continue  # Skip non-overlapping layers
        
            # Calculate properties over clipped dz
            z_vals = np.linspace(z_top_clip, z_bot_clip, npts)
            phi_vals = f_phi(z_vals)
            sigma_vals = f_sigma_v_eff(z_vals)
            delta_vals = f_delta(z_vals)
        
            phi_av = np.mean(phi_vals)
            sigma_av = np.mean(sigma_vals)
            delta_av = np.mean(delta_vals)
        
            sigma_tip = f_sigma_v_eff(z_bot_clip)
            gamma_vals = f_gamma(z_vals)
            gamma_av = np.mean(gamma_vals)
            
            Nq = np.e**(np.pi*np.tan(np.radians(phi_av)))*(np.tan(np.radians(45) + np.radians(phi_av)/2))**2
        
            # Calculate Hmax for sand
            Hmax_layer = 0.5*Nq*D*gamma_av*dz_clip**2
        
            layer_data.append({
                'z_top': z_top_clip,
                'z_bot': z_bot_clip,
                'dz': dz_clip,
                'Hmax_layer': Hmax_layer,
                'ez_layer': np.mean(z_vals)
            })
        
            # Side friction
            To = PileSurface(dz_clip, D)*delta_av*sigma_av
            Ti = PileSurface(dz_clip, D - 2*t)*delta_av*sigma_av
        
            if abs(z_bot_clip - (z0 + (L - z0))) < 1e-3:
                Tbase = np.pi/4*D**2*sigma_tip
            else:
                Tbase = 0.0
        
            Tmax = min(To + Ti, To + Tbase)
        
            # Torque induced by horizontal load
            T = Ha*rlug*np.sin(np.deg2rad(psilug))
            nhuT = T/Tmax 
            nhuV = Ha/To 
            nhuVstar = np.sqrt(nhuV**2 - nhuT**2) 
            deltastar = delta_av*(nhuVstar/nhuV) 
            
            # Constant weight
            Pile_Head = PileWeight(z0, D, t, rhows)
               
            # Vertical failure modes
            Vmax2 = Pile_Head + PileWeight(dz_clip, D, t, rhows) + PileSurface(dz_clip, D)*deltastar*sigma_av + PileSurface(dz_clip, D - 2*t)*deltastar*sigma_av
            Vmax3 = Pile_Head + PileWeight(dz_clip, D, t, rhows) + PileSurface(dz_clip, D)*deltastar*sigma_av + SoilWeight(dz_clip, D, t, gamma_av)
        
            Vmax_layer = min(Vmax2, Vmax3)
            
            # Sum vertical capacities
            Vmax_final += Vmax_layer
        
            if display > 0: print(f"Vmax_layer (sand) = {Vmax_layer:.2f} N")
            if display > 0: print(f"Vmax2 (sand) = {Vmax2:.2f} N")
            if display > 0: print(f"Vmax3 (sand) = {Vmax3:.2f} N")

    # Hmax_final and weighted ez
    for data in layer_data:
        z_top = data['z_top']
        z_bot = data['z_bot']
        Hmax_layer = data['Hmax_layer']
        ez_layer = data['ez_layer']
        dz_layer = data['dz']

        z_embedded_start = z0
        z_embedded_end = L - z0

        if z_top >= z_embedded_start and z_bot <= z_embedded_end:
            sum_ez_weighted += Hmax_layer*ez_layer
            Hmax_final += Hmax_layer
            if display > 0: print(f'Hmax_layer  = {Hmax_layer:.2f} m')

        elif z_top < z_embedded_end and z_bot > z_embedded_start:
            dz_inside = min(z_bot, z_embedded_end) - max(z_top, z_embedded_start)
            if dz_inside > 0:
                ratio = dz_inside/dz_layer
                sum_ez_weighted += Hmax_layer*ratio*ez_layer
                Hmax_final += Hmax_layer*ratio
                # print(f'ez_layer (partial) = {ez_layer:.2f} m')
   
    ez_global = sum_ez_weighted/Hmax_final 
    if display > 1: print(f'ez_global   = {ez_global:.2f} m')
    if display > 0: print(f'Hmax_final  = {Hmax_final:.2f} m')

    # Calculate coupled moment 
    M  = -Va*rlugTilt(rlug, zlug, thetalug) - Ha*(zlugTilt(rlug, zlug, thetalug) - ez_global)
    Mv = -Va*rlugTilt(rlug, zlug, thetalug)
    if display > 1: print(f"rlug_eff = {rlugTilt(rlug, zlug, thetalug):.2f} m")
    if display > 1: print(f"zlug_eff = {zlugTilt(rlug, zlug, thetalug):.2f} m")
    if display > 1: print(f"M = {M:.2f} Nm")
    
    # MH Ellipse Parameters for Clay (Kay 2014)
    # ΔφMH (piecewise based on L/D)
    if 0.5 <= lambdap < 1.125:
        delta_phi = 0.32 + 4.32*lambdap; #print(delta_phi)
    elif 1.125 <= lambdap < 2.0:
        delta_phi = 7.13 - 1.71*lambdap; #print(delta_phi)
    elif 2.0 <= lambdap <= 8.0:
        delta_phi = 2.25 - 0.25*lambdap; #print(delta_phi)
    else:
        raise ValueError('L/D out of bounds for MH ellipse.')

    phi_MH = -np.arctan(ez_global/(L - z0)) - np.deg2rad(delta_phi)
    a_MH = Np_fixed/np.cos(phi_MH)
    delta_bMH = 0.45*(lambdap)**(-0.9) if lambdap <= 1.5 else 0
    b_MH = -Np_free*np.sin(phi_MH) + delta_bMH
    if display > 1: print(f"delta_phi = {delta_phi:.2f} deg")
    if display > 1: print(f"phi_MH = {np.rad2deg(phi_MH):.2f} deg")
    if display > 1: print(f"a_MH = {a_MH:.2f}")
    if display > 1: print(f"b_MH = {b_MH:.2f}")

    # MH Ellipse Parameters for Clay (Kay 2015)
    # VH (piecewise based on L/D)
    if 0.5 <= lambdap < 1.5:
        a_VH = 9/4 + (5/3)*lambdap; 
    elif 0.5 <= lambdap < 1.25:
        b_VH = 23/4 - (13/5)*lambdap; 
    elif 1.5 <= lambdap < 20.0:     # need to set a maximum of 6 based on installation pump requirements   
        a_VH = 47/12 - (5/9)*lambdap; 
        b_VH = 50/19 - (2/19)*lambdap;
    else:
        raise ValueError('L/D ratio out of bounds for MH ellipse formulation.')
    a_VH = 0.5 + lambdap;   b_VH = 4.5 + lambdap/3
    # a_VH = 4.5 + lambdap/2; b_VH = 4.5 + lambdap/4
    if display > 0: print(f"a_VH = {a_VH:.2f}")
    if display > 0: print(f"b_VH = {b_VH:.2f}")

    # Scale VH ellipse based on vertical load ratio (Kay 2015) 
    shrink_factor = 1 - ((Va/Vmax_final)**b_VH)**(2/a_VH)

    if plot: plt.figure(figsize=(10, 5))
    theta = np.linspace(0, 2*np.pi, 400)        
    shrink_factors = np.linspace(0.0, 1.0, 5)  
    # Define colormap
    if plot: cmap = plt.colormaps['Greys']  
    norm = mcolors.Normalize(vmin=min(shrink_factors), vmax=max(shrink_factors))
    
    for s_f in shrink_factors:
        if plot: color = cmap(norm(s_f))
        x_ellipse = Hmax_final*s_f*np.cos(theta)
        y_ellipse = Vmax_final*s_f*np.sin(theta)
        H_rot = np.cos(phi_MH)*x_ellipse - np.sin(phi_MH)*y_ellipse
        M_rot = np.sin(phi_MH)*x_ellipse + np.cos(phi_MH)*y_ellipse
        if plot: plt.plot(H_rot, M_rot, color=color, alpha=0.5)

    x_ellipse_prime = Hmax_final*shrink_factor*np.cos(theta)
    y_ellipse_prime = Vmax_final*shrink_factor*np.sin(theta)
    H_rot_prime = np.cos(phi_MH)*x_ellipse_prime - np.sin(phi_MH)*y_ellipse_prime
    M_rot_prime = np.sin(phi_MH)*x_ellipse_prime + np.cos(phi_MH)*y_ellipse_prime
    Hlim = 1.2*Hmax_final  
    if plot: plt.xlim(-Hlim, Hlim)
    if plot: plt.ylim(-Hlim, Hlim)
    if plot: plt.grid(True, color='gray', linestyle='--', lw=0.5, alpha=0.8)

    # Highlight the actual one
    if plot: plt.plot(H_rot_prime, M_rot_prime, 'b', label= f'MH ellipse w/ V/Vmax = {shrink_factor:.3f}')
    if plot: plt.axhline(0, color='k', linestyle='--', lw=1.0)
    if plot: plt.axvline(0, color='k', linestyle='--', lw=1.0)

    # Plot horizontal line at constant M and Mv
    H_plot = np.linspace(min(1.3*H_rot), max(1.3*H_rot), 100)
    M_plot  = np.full_like(H_plot,  M)  # Constant moment
    Mv_plot = np.full_like(H_plot, Mv)  # Constant moment
    if plot: plt.plot(H_plot,  M_plot, 'r', lw=1.0, label='Moment line') 
    if plot: plt.plot(H_plot, Mv_plot, 'r', lw=0.5, label='Vertical moment line') 
    if plot: plt.legend(loc='lower left', fontsize='small')
   
    H_roots = horizontal_cross(H_rot_prime, M_rot_prime, M)
    Hmax_v = 0.1 
    if H_roots:
        Hmax_pos = max([r for r in H_roots if r >= 0], default=None)
        Hmax_neg = min([r for r in H_roots if r <  0], default=None)
        if M > 0 and Hmax_neg is not None:
            Hmax_v = abs(Hmax_neg)
            if plot: plt.plot(Hmax_neg, M, 'ro', label=f'Hmax,v = {Hmax_neg/1e6:.1f} MN', zorder=20)
            if plot: plt.legend(loc='lower left')
        elif M <= 0 and Hmax_pos is not None:
            Hmax_v = abs(Hmax_pos)
            if plot: plt.plot(Hmax_pos, M, 'ro', label=f'Hmax,v = {Hmax_pos/1e6:.1f} MN', zorder=20)
            if plot: plt.legend(loc='lower left')
        else:
            if display > 0: print('[WARNING] No valid Hmax crossing found for moment cut.')
    else:
        if display > 0: print('[WARNING] No intersection between moment line and ellipse.') 

    # Find relevant intercept
    H_v_roots = horizontal_cross(H_rot_prime, M_rot_prime, 0.0)
    M_v_roots = vertical_cross(H_rot_prime, M_rot_prime, 0.0)
    idx_maxH = np.argmax(H_rot_prime)
    H_at_maxH = H_rot_prime[idx_maxH]
    M_at_maxH = M_rot_prime[idx_maxH]
    idx_minM = np.argmin(M_rot_prime)
    H_at_minM = H_rot_prime[idx_minM]
    M_at_minM = M_rot_prime[idx_minM]
    
    # Plotting
    if plot:
        plt.scatter(H_v_roots[0], 0.0, s=25, facecolors='white', edgecolors='blue', 
                    marker='s',label=f'Ho ≈ {H_v_roots[0]/1e6:.1f} MN', zorder=10)
        plt.legend(loc='lower left', fontsize='small')
        plt.scatter(0.0, M_v_roots[0], s=25, facecolors='white', edgecolors='blue', 
                    marker='s', label=f'Mo ≈ {M_v_roots[0]/1e6:.1f} MNm', zorder=10)
        plt.legend(loc='lower left', fontsize='small')
        plt.scatter(H_at_maxH, M_at_maxH, s=25, facecolors='white', edgecolors='blue',
                    marker='D', label=f'Hmax ≈ {H_at_maxH/1e6:.1f} MN', zorder=10)
        plt.legend(loc='lower left', fontsize='small')
        plt.scatter(H_at_minM, M_at_minM, s=25, facecolors='white', edgecolors='blue',
                    marker='D', label=f'Mmax ≈ {M_at_minM/1e6:.1f} MNm', zorder=10)
        plt.legend(loc='lower left', fontsize='small')
    
    # Constant weight
    pile_head = PileWeight(z0, D, t, rhows)
    if display > 0: print(f"pile_head    = {pile_head:.2f} N")
    Vmax_final += pile_head
    if display > 0: print(f"Vmax_final   = {Vmax_final:.2f} N")
    
    Wp = 1.10*PileWeight(L, D, t, rhows + rhow) 
    
    # Capacity envelope
    a_VH = 0.5 + lambdap; b_VH = 4.5 + lambdap/3
    # Unity check
    UC = (Ha/Hmax_v)**a_VH + (Va/Vmax_final)**b_VH 
    if plot: plt.figure(figsize=(6, 5))
    x = np.linspace(0, 1, 100)
    y = (1 - x**b_VH)**(1/a_VH)

    # Plotting
    if plot:
        plt.figure(figsize=(6, 5))
        plt.plot(Hmax_v*x, Vmax_final*y, 'b', label='VH Envelope')
        plt.plot(Ha, Va, 'go', label='Applied load')
        plt.xlabel('Horizontal capacity (N)')
        plt.ylabel('Vertical capacity (N)')
        plt.title('VH suction pile capacity envelope')
        plt.axis([0, 1.3*max(Hmax_v, Ha), 0, 1.3*max(Vmax_final, Va)]) 
        plt.grid(True)
        plt.legend()
        plt.show()
    
    resultsSuction = {
        'Horizontal max.': Hmax_v,
        'Vertical max.': Vmax_final,
        'Unity check': UC,
        'Weight pile': Wp}

    return layers, resultsSuction

if __name__ == '__main__':

    profile_map = [
    {
        'name': 'CPT_1',
        'x': 498234, 'y': 5725141,
        'layers': [
            {
                'top': 0.0, 'bottom': 20.0,
                'soil_type': 'clay',
                'gamma_top': 8.5, 'gamma_bot': 8.5,
                'Su_top': 2.4, 'Su_bot': 30.3}]
            # {
            #     'top': 10.0, 'bottom': 20.0,
            #     'soil_type': 'clay',
            #     'gamma_top': 8.5, 'gamma_bot': 8.5,
            #     'Su_top': 13.95, 'Su_bot': 30.3}]
            # {
            #     'top': 30.0, 'bottom': 36.0,
            #     'soil_type': 'clay',
            #     'gamma_top': 9.0, 'gamma_bot': 9.5,
            #     'Su_top': 75, 'Su_bot': 100},
            # {
            #     'top': 36.0, 'bottom': 55.0,
            #     'soil_type': 'clay',
            #     'gamma_top': 9.5, 'gamma_bot': 9.5,
            #     'Su_top': 100, 'Su_bot': 100}]
        }
    ]


    # Pile and load properties
    D = 3.34                          # Pile diameter (m)
    L = 20.0                          # Pile length (m)
    zlug = (2/3)*L                    # Lug depth (m)
    theta = 5                         # Tilt misalignment angle (deg)
    psi = 7.5                         # Twist misalignment angle (deg)
    Ha = 1e6                          # Applied horizontal load (N)
    Va = 5.7e6                        # Applied vertical load (N)

    # Calculate
    layers, resultsSuction = getCapacitySuction(
        profile_map, 'CPT_1',         # Soil properties and location of the pile
        D, L, zlug,                   # Pile geometrical properties
        Ha, Va,                       # Pile loading conditions   
        thetalug=theta, psilug=psi,   # Pile misaligment tolerances 
        plot=True,
        display=1
    )

    # print('\n--- Suction Pile Capacity Results ---')
    # print(f"Hmax_final = {resultsSuction['Hmax_final']:.2f} N")
    # print(f"Vmax_final = {resultsSuction['Vmax_final']:.2f} N")
    # print(f"Unity check (UC) = {resultsSuction['UnityCheck']:.4f}")
    # print(f"Total Moment (M_total) = {resultsSuction['M_total']:.2f} Nm")

    # plot_suction(layers, L, D, z0 = layers[0]['top'], zlug=zlug)
