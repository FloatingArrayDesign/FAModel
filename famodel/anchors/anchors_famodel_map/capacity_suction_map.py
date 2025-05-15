
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from .capacity_soils_map import clay_profile, sand_profile
from .capacity_plots_map import plot_suction


def PileSurface(Len, Dia):
    return np.pi*Dia*Len

def PileWeight(Len, Dia, tw, rho):
    return ((np.pi/4)*((Dia**2 - (Dia - 2*tw)**2)*Len + (np.pi/4)*Dia**2*tw))*rho

def SoilWeight(Len, Dia, tw, gamma_soil):
    return (np.pi/4)*(Dia - 2*tw)**2*Len*gamma_soil

def rlugTilt(r, z, theta):
    return r*np.cos(np.deg2rad(theta)) - z*np.sin(np.deg2rad(theta))

def zlugTilt(r, z, theta):
    return r*np.sin(np.deg2rad(theta)) + z*np.cos(np.deg2rad(theta))

def getCapacitySuction(profile_map, location_name, D, L, zlug, Ha, Va, thetalug=5, psilug=7.5, plot=True):
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
    
    # Retrieve soil layers from map
    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']
    
    z0 = layers[0]['top']            # Mudline elevation    
    lambdap = (L - z0)/D             # Suction pile slenderness ratio
    t = (6.35 + D*20)/1e3            # Suction pile wall thickness (m), API RP2A-WSD
    rlug = D/2                       # Radial position of the lug
    thetalug = 5                     # Angle of tilt misaligment, default is 5. (deg)
    psilug = 7.5                     # Angle of twist misaligment, default is 7.5. (deg)
    rhows = 66.90e3                  # Submerged steel specific weight (N/m3)
    rhow = 10e3                      # Water specific weight (N/m3) 
    
    Np_fixed = 10.25
    Np_free = 4.0
    Nc = min(6.2*(1 + 0.34*np.arctan(lambdap)), 9)
    
    # Initialize
    sum_ez_weighted = 0.0
    sum_Hmax = 0.0
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
            dz_clip = z_bot_clip - z_top_clip; # print(f'dz_clip  = {dz_clip:.2f} m')
    
            if dz_clip <= 0:
                continue  # Skip layers fully above or below
    
            # Calculate properties over clipped dz
            z_vals = np.linspace(z_top_clip, z_bot_clip, npts)
            Su_vals = f_Su(z_vals)
            Su_total = np.trapz(Su_vals, z_vals)
            Su_moment = np.trapz(Su_vals*z_vals, z_vals)
    
            Su_av_z = Su_total/dz_clip; # print(f'Su_av_z  = {Su_av_z:.2f} Pa')
            ez_layer = Su_moment/Su_total; 
            Su_bot = f_Su(z_bot_clip)
            gamma_vals = f_gamma(z_vals)
            gamma_av = np.mean(gamma_vals)
    
            # Calculate Hmax for clay
            Hmax_layer = Np_fixed*D*dz_clip*Su_av_z; print(f'Su_av_z  = {Su_av_z:.2f} Pa')
    
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
            alphastar = alpha_av*(nhuVstar/nhuV); print(f"alphastar   = {alphastar:.3f}")
            
            # Constant weight
            Pile_Head = PileWeight(z0, D, t, rhows)
            
            # Vertical failure modes
            if np.isclose(z_bot_clip, z0 + (L - z0), atol=0.1):
                Vmax1 = PileWeight(dz_clip, D, t, rhows) + PileSurface(dz_clip, D)*alphastar*Su_av_z + Nc*Su_bot*(np.pi/4)*D**2
            else:
                Vmax1 = np.inf  # No tip resistance unless at tip
            
            Vmax2 = PileWeight(dz_clip, D, t, rhows) + PileSurface(dz_clip, D)*alphastar*Su_av_z + PileSurface(dz_clip, D - 2*t)*alphastar*Su_av_z
            Vmax3 = PileWeight(dz_clip, D, t, rhows) + PileSurface(dz_clip, D)*alphastar*Su_av_z + SoilWeight(dz_clip, D, t, gamma_av)
            
            Vmax_layer = min(Vmax1, Vmax2, Vmax3)
            
            # Sum vertical capacities
            Vmax_final += Vmax_layer
            
            # Print layer debug info
            print(f"Vmax_layer    = {Vmax_layer:.2f} N")
            print(f"Vmax1   = {Vmax1:.2f} N")
            print(f"Vmax2   = {Vmax2:.2f} N")
            print(f"Vmax3   = {Vmax3:.2f} N")
            
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
               
            # Vertical failure modes
            Vmax2 = Pile_Head + PileWeight(dz_clip, D, t, rhows) + PileSurface(dz_clip, D)*deltastar*sigma_av + PileSurface(dz_clip, D - 2*t)*deltastar*sigma_av
            Vmax3 = Pile_Head + PileWeight(dz_clip, D, t, rhows) + PileSurface(dz_clip, D)*deltastar*sigma_av + SoilWeight(dz_clip, D, t, gamma_av)
        
            Vmax_layer = min(Vmax2, Vmax3)
            
            # Sum vertical capacities
            Vmax_final += Vmax_layer
        
            print(f"Vmax_layer (sand) = {Vmax_layer:.2f} N")
            print(f"Vmax2 (sand) = {Vmax2:.2f} N")
            print(f"Vmax3 (sand) = {Vmax3:.2f} N")

    # sum Hmax and weighted ez
    for data in layer_data:
        z_top = data['z_top']
        z_bot = data['z_bot']
        Hmax_layer = data['Hmax_layer']
        ez_layer = data['ez_layer']
        dz_layer = data['dz']

        z_embedded_start = z0
        z_embedded_end = z0 + (L - z0)

        if z_top >= z_embedded_start and z_bot <= z_embedded_end:
            sum_ez_weighted += Hmax_layer*ez_layer
            sum_Hmax += Hmax_layer
            # print(f'ez_layer (full) = {ez_layer:.2f} m')

        elif z_top < z_embedded_end and z_bot > z_embedded_start:
            dz_inside = min(z_bot, z_embedded_end) - max(z_top, z_embedded_start)
            if dz_inside > 0:
                ratio = dz_inside/dz_layer
                sum_ez_weighted += Hmax_layer*ratio*ez_layer
                sum_Hmax += Hmax_layer * ratio
                # print(f'ez_layer (partial) = {ez_layer:.2f} m')

    ez_global = sum_ez_weighted/sum_Hmax 
    # print(f'sum_ez_weighted = {sum_ez_weighted:.2f}')
    print(f'ez_global = {ez_global:.2f} m')
    print(f'Hmax  = {sum_Hmax:.2f} m')

    # Calculate moment and Hmax_final
    M_total = -Va*rlugTilt(rlug, zlug, thetalug) - Ha*(zlugTilt(rlug, zlug, thetalug) - ez_global)
    # print(f"rlug_eff = {rlugTilt(rlug, zlug, thetalug):.2f} m")
    # print(f"zlug_eff = {zlugTilt(rlug, zlug, thetalug):.2f} m")
    print(f"M_total = {M_total:.2f} Nm")
    
    # ΔφMH from Kay 2014
    if 0.5 <= lambdap < 1.125:
        delta_phi = 0.32 + 4.32*lambdap; #print(delta_phi)
    elif 1.125 <= lambdap < 2.0:
        delta_phi = 7.13 - 1.71*lambdap; #print(delta_phi)
    elif 2.0 <= lambdap <= 6.0:
        delta_phi = 4.55 - 0.425*lambdap; #print(delta_phi)
    else:
        raise ValueError('L/D out of bounds for MH ellipse.')

    phi_MH = -np.arctan(ez_global/(L - z0)) - np.deg2rad(delta_phi)
    a_MH = Np_fixed/np.cos(phi_MH)
    delta_bMH = 0.45*(lambdap)**(-0.9) if lambdap <= 1.5 else 0
    b_MH = -Np_free*np.sin(phi_MH) + delta_bMH
    print('M cos(phi)/a_MH =', (M_total*np.cos(phi_MH))/a_MH)
    print('M sin(phi)/b_MH =', (M_total*np.sin(phi_MH))/b_MH)

    def f(H_var):
        term1 = ((M_total*np.cos(phi_MH) + H_var*np.sin(phi_MH))/a_MH)**2
        term2 = ((M_total*np.sin(phi_MH) - H_var*np.cos(phi_MH))/b_MH)**2
        return term1 + term2 - 1

    try:
        Hmax_final = max(fsolve(f, sum_Hmax*0.8)[0], 0.0)
    except:
        Hmax_final = 0.0

    print(f"Hmax_final (MH ellipse) = {Hmax_final:.2f} N")
    
    # Constant weight
    pile_head = PileWeight(z0, D, t, rhows); print(f"pile_head    = {pile_head:.2f} N")
    Vmax_final += pile_head; print(f"Vmax_final    = {Vmax_final:.2f} N")
    
    # Unity check
    UC = (Ha/Hmax_final)**(0.5 + lambdap) + (Va/Vmax_final)**(4.5 + lambdap/3) if Hmax_final and sum_Hmax else np.inf

    # Plotting
    if plot:
        x = np.linspace(0, 1, 100)
        y = (1 - x**(4.5 + lambdap/3))**(1/(0.5 + lambdap))
        
        plt.figure(figsize=(6, 5))
        plt.plot(Hmax_final*x, Vmax_final*y, 'b', label='VH Envelope')
        plt.plot(Ha, Va, 'ro', label='Applied Load')
        plt.xlabel('Horizontal Capacity (N)')
        plt.ylabel('Vertical Capacity (N)')
        plt.title('VH suction pile capacity envelope')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

    resultsSuction = {
        'Horizontal max.': Hmax_final,
        'Vertical max.': Vmax_final,
        'Unity check': UC,
        # 'Weight pile': Wp,
        # 'Weight soil': Wsoil,
        't': t
    }

    return layers, resultsSuction

if __name__ == '__main__':

    profile_map = [
    {
        'name': 'CPT_1',
        'x': 498234, 'y': 5725141,
        'layers': [
            {
                'top': 2.0, 'bottom': 4.0,
                'soil_type': 'clay',
                'gamma_top': 8.0, 'gamma_bot': 8.5,
                'Su_top': 25, 'Su_bot': 50},
            {
                'top': 4.0, 'bottom': 8.0,
                'soil_type': 'clay',
                'gamma_top': 8.5, 'gamma_bot': 9.0,
                'Su_top': 50, 'Su_bot': 75},
            {
                'top': 8.0, 'bottom': 16.0,
                'soil_type': 'clay',
                'gamma_top': 9.0, 'gamma_bot': 9.5,
                'Su_top': 75, 'Su_bot': 100},
            {
                'top': 16.0, 'bottom': 25.0,
                'soil_type': 'clay',
                'gamma_top': 9.5, 'gamma_bot': 9.5,
                'Su_top': 100, 'Su_bot': 100}]
        }
    ]


    # Pile and load properties
    D = 2.5                           # Pile diameter (m)
    L = 10.0                          # Pile length (m)
    zlug = 8.0                        # Lug depth (m)
    theta = 5                         # Tilt misalignment angle (deg)
    psi = 7.5                         # Twist misalignment angle (deg)
    Ha = 6e6                          # Applied horizontal load (N)
    Va = 2e6                          # Applied vertical load (N)

    # Calculate
    layers, resultsSuction = getCapacitySuction(
        profile_map, 'CPT_1',         # Soil properties and location of the pile
        D, L, zlug,                   # Pile geometrical properties
        Ha, Va,                       # Pile loading conditions   
        thetalug=theta, psilug=psi,   # Pile misaligment tolerances 
        plot=True
    )

    # print('\n--- Suction Pile Capacity Results ---')
    # print(f"Hmax_final = {resultsSuction['Hmax_final']:.2f} N")
    # print(f"Vmax_final = {resultsSuction['Vmax_final']:.2f} N")
    # print(f"Unity check (UC) = {resultsSuction['UnityCheck']:.4f}")
    # print(f"Total Moment (M_total) = {resultsSuction['M_total']:.2f} Nm")

    plot_suction(layers, L, D, z0 = layers[0]['top'], zlug=zlug, title='Suction Pile and Soil Layers')
