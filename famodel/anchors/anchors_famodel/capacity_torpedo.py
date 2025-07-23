
import numpy as np
import matplotlib.pyplot as plt
from .support_soils import clay_profile
from .support_plots import plot_torpedo

def getCapacityTorpedo(profile_map, location_name, D1, D2, L1, L2, zlug, ballast, Ha, Va, plot=False):
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

    # Retrieve soil layers from map
    profile_entry = next(p for p in profile_map if p['name'] == location_name)
    layers = profile_entry['layers']

    L = L1 + L2
    t = (6.35 + D2*20)/1e3
    rhows = 66.90e3
    rhow = 10e3

    def PileWeight(Len1, Len2, Dia1, Dia2, tw, rho):
        return ((np.pi/4)*(Dia1**2 - (Dia1 - 2*tw)**2)*(Len1 + Len2) + 4*Len2*Dia2*tw)*rho

    def PileWingedSurface(length, diameter):
        return np.pi*diameter*length

    def PileShaftSurface(length, diameter1, diameter2):
        return 8*length*(diameter1 - diameter2)

    z_start = zlug
    z_wing = zlug + L1
    z_end = zlug + L

    layer_data = []
    Vmax_total = 0.0
    
    # Profile check points
    npts = 10

    for layer in layers:
        if layer['soil_type'] != 'clay':
            raise ValueError('Torpedo pile capacity model only supports clay soils.')

        z_layer_top = layer['top']
        z_layer_bot = layer['bottom']

        z_clip_top = max(z_layer_top, z_start)
        z_clip_bot = min(z_layer_bot, z_end)

        if z_clip_bot <= z_clip_top:
            continue

        segments = []
        if z_clip_bot <= z_wing:
            segments.append((z_clip_top, z_clip_bot, D1))
        elif z_clip_top >= z_wing:
            segments.append((z_clip_top, z_clip_bot, D2))
        else:
            segments.append((z_clip_top, z_wing, D1))
            segments.append((z_wing, z_clip_bot, D2))

        for z_seg_top, z_seg_bot, D in segments:
            dz_seg = z_seg_bot - z_seg_top
            if dz_seg <= 0:
                continue

            profile = [
                [z_seg_top, layer['Su_top'], layer['gamma_top']],
                [z_seg_bot, layer['Su_bot'], layer['gamma_bot']]
            ]
            z_ref, f_Su, _, f_gamma, f_alpha = clay_profile(profile)

            z_vals = np.linspace(z_seg_top, z_seg_bot, npts)
            Su_vals = f_Su(z_vals)
            alpha_vals = np.array([f_alpha(z) for z in z_vals])

            Su_total = np.trapz(Su_vals, z_vals)
            Su_moment = np.trapz(z_vals*Su_vals, z_vals)
            print("xxxxxxxxxxxxxxxxxxxxxxxxx")
            Su_av_z = Su_total/dz_seg
            print(f"Su_av_z = {Su_av_z:.2f} Pa")
            ez_layer = Su_moment /Su_total
            print(f"dz_seg = {dz_seg:.2f} m")
            print(f"ez_layer = {ez_layer:.2f} m")
            alpha_av = np.mean(alpha_vals)
            print(f"alpha_av = {alpha_av:.2f}")

            Np_free = 3.45
            Hmax_layer = Np_free*dz_seg*D*Su_av_z
            print(f"Hmax_layer = {Hmax_layer:.2f} N")
            print(f"D = {D:.2f} m")

            surface_area = PileWingedSurface(dz_seg, D) if D == D1 else PileShaftSurface(dz_seg, D1, D2)
            Vmax_layer = surface_area*alpha_av*Su_av_z
            Vmax_total += Vmax_layer
            print(f"Vmax_layer = {Vmax_layer:.2f} N")
            
            layer_data.append({
                'z_top': z_seg_top,
                'z_bot': z_seg_bot,
                'dz': dz_seg,
                'Hmax_layer': Hmax_layer,
                'ez_layer': ez_layer,
                'Su_av_z': Su_av_z,
                'D_used': D
            })

    if not layer_data:
        raise ValueError('No overlapping clay layers within pile depth.')

    sum_Hmax = 0.0
    sum_ez_weighted = 0.0

    for data in layer_data:
        z_top = data['z_top']
        z_bot = data['z_bot']
        Hmax_layer = data['Hmax_layer']
        ez_layer = data['ez_layer']
        dz_layer = data['dz']

        z_embedded_start = zlug
        z_embedded_end = zlug + L

        if z_top >= z_embedded_start and z_bot <= z_embedded_end:
            sum_ez_weighted += Hmax_layer*ez_layer
            sum_Hmax += Hmax_layer
        elif z_top < z_embedded_end and z_bot > z_embedded_start:
            dz_inside = min(z_bot, z_embedded_end) - max(z_top, z_embedded_start)
            if dz_inside > 0:
                ratio = dz_inside/dz_layer
                sum_ez_weighted += Hmax_layer*ratio*ez_layer
                sum_Hmax += Hmax_layer * ratio

    ez_global = sum_ez_weighted/sum_Hmax
    print(f'ez_global = {ez_global:.2f} m')    
    print(f'sum_Hmax = {sum_Hmax:.2f} N')

    Vmax_total += PileWeight(L1, L2, D1, D2, t, rhows) + ballast
    Wp = 1.10 * PileWeight(L1, L2, D1, D2, t, rhows + rhow) + ballast

    ez_ratio = (ez_global - zlug)/L; print(f"ez_ratio = {ez_ratio:.2f} m")
    
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

    # Assign aVH and bVH based on ez_su/L
    if np.isclose(ez_ratio, 2/3, atol=0.05):
        aVH = 0.5 + L/Dstar
        bVH = 4.5 - L/(3*Dstar)
        mode = 'deep mobilization (2/3)'
    elif 0.40 <= ez_ratio <= 0.75:
        aVH = 4.5 + L/(2*Dstar)
        bVH = 3.5 - L/(4*Dstar)
        mode = 'moderate mobilization (1/2 – 3/4)'
    # else:
    #     aVH = 4.0
    #     bVH = 4.0
    #     mode = 'default exponents (fallback)'
    print(f'Interaction exponents set to aVH = {aVH:.2f}, bVH = {bVH:.2f} [{mode}]')
    
    UC = (Ha/sum_Hmax)**aVH + (Va/Vmax_total)**bVH

    if plot:
        deg = np.linspace(0, 90, 20)
        x = np.cos(np.deg2rad(deg))
        y = (1 - x**bVH)**(1/aVH)
        X = sum_Hmax*x
        Y = Vmax_total*y

        plt.figure(figsize=(6, 5))
        plt.plot(X, Y, color='blue', label='VH Envelope')
        plt.plot(Ha, Va, 'o', color='red', label='Load Point')
        plt.xlabel('Horizontal Capacity (N)')
        plt.ylabel('Vertical Capacity (N)')
        plt.title('VH torpedo pile capacity envelope')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

    resultsTorpedo = {
        'Horizontal max.': sum_Hmax,
        'Vertical max.': Vmax_total,
        'Unity check': UC,
        'Weight pile': Wp,
        'ez_global': ez_global,
        'layer_data': layer_data}

    return layers, resultsTorpedo

if __name__ == '__main__':

    profile_map = [
        {
            'name': 'CPT_1',
            'x': 498234, 'y': 5725141,
            'layers': [
                {
                    'top': 0.0, 'bottom': 20.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.0, 'gamma_bot': 8.5,
                    'Su_top': 50, 'Su_bot': 70},
                {
                    'top': 20.0, 'bottom': 25.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.5, 'gamma_bot': 8.5,
                    'Su_top': 80, 'Su_bot': 100},
                {
                    'top': 25.0, 'bottom': 50.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.5, 'gamma_bot': 9.0,
                    'Su_top': 125, 'Su_bot': 150}]
        }
    ]

    D1 = 3.0
    D2 = 1.5
    L1 = 11.0
    L2 = 5.0
    zlug = 15.0
    ballast = 10000
    Ha = 6.0e6
    Va = 8.0e6

    layers, results = getCapacityTorpedo(profile_map, 'CPT_1', D1, D2, L1, L2, zlug, ballast, Ha, Va)

    # print("\n--- Torpedo Pile Capacity Results ---")
    # for key, val in results.items():
    #     if key != 'layer_data':
    #         print(f"{key}: {val:.2f}")
            
    plot_torpedo(layers, D1, D2, L1, L2, z0 = layers[0]['top'], zlug=zlug)
