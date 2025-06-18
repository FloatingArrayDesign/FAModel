
import numpy as np
import matplotlib.pyplot as plt

def plot_pile(layers, y, z, D, L, z0=None, zlug=None, hinge_location=None):
    fig, ax = plt.subplots(figsize=(5, 5))

    lambdap = L / D
    if lambdap >= 4:
        xmax = 5*D
    elif lambdap <= 2:
        xmax = 2*D
    else:
        xmax = 3*D

    # Normalize color ranges based on max values
    max_Su  = max((layer.get( 'Su_top', 0) for layer in layers if layer['soil_type'] == 'clay'), default=1)
    max_phi = max((layer.get('phi_top', 0) for layer in layers if layer['soil_type'] == 'sand'), default=1)
    max_UCS = max((layer.get('UCS_top', 0) for layer in layers if layer['soil_type'] in ['rock', 'weak_rock']), default=1)

    seen_labels = set()
    for layer in layers:
        z_top = layer['top']
        z_bot = layer['bottom']
        soil = layer['soil_type']

        if soil == 'clay':
            Su = layer.get('Su_top', 25)
            gamma = layer.get('gamma_top', 10)
            color_fill = plt.cm.Oranges(Su/max_Su)
            label_soil = f'Su = {Su:.0f} kPa, γ = {gamma:.1f} kN/m³'
        elif soil == 'sand':
            phi = layer.get('phi_top', 30)
            gamma = layer.get('gamma_top', 10)
            color_fill = plt.cm.YlOrBr(phi/max_phi)
            label_soil = f'ϕ = {phi:.0f}°, γ = {gamma:.1f} kN/m³'
        elif soil in ['rock', 'weak_rock']:
            UCS = layer.get('UCS_top', 1.0)
            Em = layer.get('Em_top', 50)
            color_fill = plt.cm.Greys(UCS/max_UCS)
            label_soil = f'UCS = {UCS:.2f} MPa, Em = {Em:.1f} MPa'
        else:
            color_fill = 'gray'
            label_soil = soil.capitalize()

        if label_soil not in seen_labels:
            ax.axhspan(z_bot, z_top, color=color_fill, alpha=0.4, label=label_soil)
            seen_labels.add(label_soil)
        else:
            ax.axhspan(z_bot, z_top, color=color_fill, alpha=0.4)
    
    # Draw pile geometry
    pile = plt.Rectangle((-D/2, 0), D, L, edgecolor='k', facecolor='none', lw=2, label='Driven pile')
    ax.add_patch(pile)
    
    # Padeye marker
    if zlug is not None:
        ax.plot(D/2, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')
    
    if hinge_location is not None and hinge_location != -1:
        ax.plot(y[hinge_location], z[hinge_location], 'o', color='black', label='Plastic hinge')
        
    # Reference lines
    ax.axhline(z0, color='b', linestyle='--', lw=1.5, label='Mudline')
    ax.axhline( L, color='r', linestyle='--', lw=1.5, label='Pile tip')

    ax.set_xlabel('Lateral displacement (m)')
    ax.set_ylabel('Depth (m)')
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim([L + 5, -2])
    ax.grid(ls='--')
    ax.legend()
    ax.set_title('Pile Deflection Profile')
    plt.tight_layout()
    plt.show()
    
def plot_suction(layers, L, D, z0=None, zlug=None, title='Suction Pile and Soil Layers'):
    '''Plot the soil profile and a suction pile geometry using updated profile_map structure.

    Parameters:
    ----------
        layers : list of dicts 
            Each dict has 'top', 'bottom', 'soil_type', and top-of-layer properties
        L : float 
            Embedded length (m)
        D : float 
            Pile diameter (m)
        zlug : float
            Padeye depth (m, referenced to pile head = 0)
        title : string 
            Plot title
    '''
    fig, ax = plt.subplots(figsize=(8, 5))
    xmax = 2*D
    
    # Normalize for each soil type
    max_Su  = max((layer.get( 'Su_top', 0) for layer in layers if layer['soil_type'] == 'clay'), default=1)
    max_phi = max((layer.get('phi_top', 0) for layer in layers if layer['soil_type'] == 'sand'), default=1)

    seen_labels = set()
    for layer in layers:
        z_top = layer['top']
        z_bot = layer['bottom']
        soil = layer['soil_type']

        if soil == 'clay':
            Su = layer.get('Su_top', 25)
            gamma = layer.get('gamma_top', 10)
            color_fill = plt.cm.Oranges(Su / max_Su)
            label_soil = f'Su = {Su:.0f} kPa, γ = {gamma:.1f} kN/m³'
        elif soil == 'sand':
            phi = layer.get('phi_top', 30)
            gamma = layer.get('gamma_top', 10)
            color_fill = plt.cm.YlOrBr(phi / max_phi)
            label_soil = f'ϕ = {phi:.0f}°, γ = {gamma:.1f} kN/m³'
        else:
            color_fill = 'gray'
            label = soil.capitalize()

        if label_soil not in seen_labels:
            ax.axhspan(z_bot, z_top, color=color_fill, alpha=0.4, label=label_soil)
            seen_labels.add(label_soil)
        else:
            ax.axhspan(z_bot, z_top, color=color_fill, alpha=0.4)

    # Draw pile geometry
    x_left = -D/2; x_right = D/2
    z_top = 0; z_bot = L
    ax.plot([  x_left, x_left], [z_top, z_bot], color='k', lw=2.0, label='Suction pile')
    ax.plot([x_right, x_right], [z_top, z_bot], color='k', lw=2.0)
    ax.plot([ x_left, x_right], [z_top, z_top], color='k', lw=2.0)
    
    # Padeye marker
    if zlug is not None:
        ax.plot(x_right, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')

    # Reference lines
    ax.axhline(z0, color='b', linestyle='--', lw=1.5, label='Mudline')
    ax.axhline( L, color='r', linestyle='--', lw=1.5, label='Pile tip')

    ax.set_xlabel('Horizontal extent (m)')
    ax.set_ylabel('Depth (m)')
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(L + 2*D, -D)
    ax.set_title(title)
    ax.grid()
    ax.legend()
    plt.tight_layout()
    plt.show()

def plot_torpedo(layers, D1, D2, L1, L2, z0, zlug, title='Torpedo Pile and Soil Layers'):
    '''Plot the soil layers and geometry of a torpedo pile using absolute depth for soil and pile head at z=0.

    Parameters:
    ----------
        profile : list of tuples  
            Each tuple is (z, Su, gamma)
        soil_type : string
            Select soil condition, 'clay' 
        D1 : float 
            Wing diameter (at tip) (m)
        D2 : float
            Shaft diameter (at head)
        L1 : float
            Winged length (m)
        L2 : float
            Shaft length (m)
        title : str
            Plot title
    '''
    fig, ax = plt.subplots(figsize=(7, 7))

    xmax = 5*max(D1, D2)
    z1 = zlug + L1            # Wing-shaft interface
    z_tip = zlug + L1 + L2    # Total depth of pile

    # Normalize color scale for Su
    max_Su = max((layer.get('Su_top', 0) for layer in layers if layer['soil_type'] == 'clay'), default=1)

    seen_labels = set()
    for layer in layers:
        z_top = layer['top']
        z_bot = layer['bottom']
        soil = layer['soil_type']

        if soil == 'clay':
            Su = layer.get('Su_top', 25)
            gamma = layer.get('gamma_top', 10)
            color_fill = plt.cm.Oranges(Su/max_Su)
            label_soil = f'Su = {Su:.0f} kPa, γ = {gamma:.1f} kN/m³'
        else:
            color_fill = 'gray'
            label_soil = soil.capitalize()

        if label_soil not in seen_labels:
            ax.axhspan(z_bot, z_top, color=color_fill, alpha=0.4, label=label_soil)
            seen_labels.add(label_soil)
        else:
            ax.axhspan(z_bot, z_top, color=color_fill, alpha=0.4)

    Df = 0.5*(D1 - D2)

    # Draw lateral fins
    ax.add_patch(plt.Rectangle((-D2/2 - Df, zlug), Df, L1, edgecolor='k', facecolor='none', lw=2, label='Fins'))
    ax.add_patch(plt.Rectangle((      D2/2, zlug), Df, L1, edgecolor='k', facecolor='none', lw=2))

    # Vertical fin
    ax.add_patch(plt.Rectangle((-0.05*Df, zlug), 0.1*Df, L1, edgecolor='k', lw=2))

    # Shaft section
    ax.add_patch(plt.Rectangle((-D2/2, zlug), D2, L1 + L2, edgecolor='k', facecolor='none', lw=2, label='Shaft section'))

    # Padeye
    if zlug is not None:
        ax.plot(0, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')
        
    # Reference lines
    ax.axhline(   z0, color='b', linestyle='--', lw=1.5, label='Mudline')
    ax.axhline(z_tip, color='r', linestyle='--', lw=1.5, label='Pile tip')

    # Axis and formatting
    zmin = min(zlug, z0) - 3*D2
    zmax = max(z_tip, z0) + 3*D2

    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(zmax, zmin)
    ax.set_xlabel('Horizontal extent (m)')
    ax.set_ylabel('Depth (m)')
    ax.set_title(title)
    ax.grid()
    ax.legend()
    plt.tight_layout()
    plt.show()
 
def plot_helical(layers, D, L, d, z0, zlug, n_helix=1, spacing=1.0, title='Helical Pile and Soil Layers'):
    '''Plot a helical pile in layered soil with shaft and angled helices, starting at zlug.

    Parameters:
    ----------
    layers : list of dicts
        Each layer has 'top', 'bottom', 'soil_type', and top-of-layer properties.
    D : float
        Helix diameter (m)
    L : float
        Shaft length (m)
    d : float
        Shaft diameter (m)
    zlug : float
        Depth of pile head
    n_helix : int
        Number of helices (typically 1)
    spacing : float
        Vertical spacing between helices (m)
    title : str
        Plot title
    '''
    fig, ax = plt.subplots(figsize=(5, 6))

    lambdap = L/D
    if lambdap >= 4:
        xmax = 5*D
    elif lambdap <= 2:
        xmax = 2*D
    else:
        xmax = 3*D

    z_tip = zlug + L

    # Normalize color scales
    max_Su  = max((layer.get('Su_top', 0) for layer in layers if layer['soil_type'] == 'clay'), default=1)
    max_phi = max((layer.get('phi_top', 0) for layer in layers if layer['soil_type'] == 'sand'), default=1)

    seen_labels = set()
    for layer in layers:
        z_top = layer['top']
        z_bot = layer['bottom']
        soil = layer['soil_type']

        if soil == 'clay':
            Su = layer.get('Su_top', 25)
            gamma = layer.get('gamma_top', 10)
            color_fill = plt.cm.Oranges(Su/max_Su)
            label_soil = f'Su = {Su:.0f} kPa, γ = {gamma:.1f} kN/m³'
        elif soil == 'sand':
            phi = layer.get('phi_top', 30)
            gamma = layer.get('gamma_top', 10)
            color_fill = plt.cm.YlOrBr(phi/max_phi)
            label_soil = f'ϕ = {phi:.0f}°, γ = {gamma:.1f} kN/m³'
        else:
            color = 'gray'
            label_soil = soil.capitalize()

        if label_soil not in seen_labels:
            ax.axhspan(z_bot, z_top, color=color_fill, alpha=0.4, label=label_soil)
            seen_labels.add(label_soil)
        else:
            ax.axhspan(z_bot, z_top, color=color_fill, alpha=0.4)

    # Shaft as vertical rectangle
    ax.add_patch(plt.Rectangle((-d/2, 0), d, L, edgecolor='k', facecolor='none', lw=2, label='Shaft'))

    # Draw angled helices at base (start from zlug + L - D)
    base_helix = L - D
    helix_height = 0.5*D
    for i in range(n_helix):
        z_helix = base_helix - i*spacing
        if z_helix < zlug:
            break
        x_helix = [-D/2, D/2]
        y_helix = [z_helix - helix_height/2, z_helix + helix_height/2]
        ax.plot(x_helix, y_helix, color='k', lw=2, label='Helix' if i == 0 else None)

    # Padeye marker on right shaft
    if zlug is not None:
        ax.plot(d/2, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')

    # Reference lines
    ax.axhline( z0, color='b', linestyle='--', lw=1.5, label='Mudline')
    ax.axhline(  L, color='k', linestyle='--', lw=1.5, label='Pile tip')

    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(L + D, min(zlug - D, min(layer['top'] for layer in layers) - 2))
    ax.set_xlabel('Horizontal extent (m)')
    ax.set_ylabel('Depth (m)')
    ax.set_title(title)
    ax.grid()
    ax.legend()
    plt.tight_layout()
    plt.show()

def plot_plate(layers, B, L, z0, zlug, beta, title='Plate Anchor in Layered Soil'):
    '''Plot soil layers and an inclined plate anchor centered at zlug.

    Parameters:
    ----------
        layers : list of dicts 
            Each dict has 'top', 'bottom', 'soil_type', and top-of-layer properties
        B : float
            Plate width (m)
        L : float
            Plate length (m) 
        z0 : float
            Mudline depth (m)
        zlug : float
            Center embedment of the plate (m)
        beta : float
            Inclination angle of plate (deg)
        title : str
            Plot title
    '''
    fig, ax = plt.subplots(figsize=(5, 5))
    xmax = 3*B

    # Normalize color scales
    max_Su  = max((layer.get('Su_top', 0) for layer in layers if layer['soil_type'] == 'clay'), default=1)

    seen_labels = set()
    for layer in layers:
        z_top = layer['top']
        z_bot = layer['bottom']
        soil = layer['soil_type']

        if soil == 'clay':
            Su = layer.get('Su_top', 25)
            color_fill = plt.cm.Oranges(Su / max_Su)
            label_soil = f'Su = {Su:.0f} kPa'
        else:
            color_fill = 'gray'
            label_soil = soil.capitalize()

        if label_soil not in seen_labels:
            ax.axhspan(z_bot, z_top, color=color_fill, alpha=0.4, label=label_soil)
            seen_labels.add(label_soil)
        else:
            ax.axhspan(z_bot, z_top, color=color_fill, alpha=0.4)

    # Inclined plate geometry
    dx = (B/2)*np.cos(np.deg2rad(beta))
    dz = (B/2)*np.sin(np.deg2rad(beta))
    plate_x = [-dx, dx]
    plate_z = [zlug - dz, zlug + dz]
    ax.plot(plate_x, plate_z, color='k', lw=1.5, label='Inclined plate')
       
    # Padeye marker
    if zlug is not None:
        ax.plot(0, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')

    # Reference lines
    ax.axhline(        z0, color='b', linestyle='--', lw=1.5, label='Mudline')
    ax.axhline(plate_z[1], color='r', linestyle='--', lw=1.5, label='Plate tip')

    # View limits
    zmin = min(zlug - dz, z0) - 2*B
    zmax = max(zlug + dz, z0) + 2*B

    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(zmax, zmin)
    ax.set_xlabel("Horizontal extent (m)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(title)
    ax.legend(loc='lower right')
    ax.grid(True)
    plt.tight_layout()
    plt.show()

    
def plot_load(layers, drag_values, depth_values, Tm, thetam, Ta, thetaa, zlug):
    '''Plot the inverse catenary profile of a mooring line with layered soil and applied load vectors.

    Parameters
    ----------
    layers : list of dicts
        Soil layers with 'top', 'bottom', 'soil_type', and property fields
    drag_values : list or array
        Horizontal drag distances of the mooring line (m)
    depth_values : list or array
        Vertical depths of the mooring line (m)
    Tm : float
        Magnitude of load applied at the mudline (N)
    thetam : float
        Inclination angle of mudline load (deg)
    Ta : float
        Magnitude of load applied at the padeye (N)
    thetaa : float
        Inclination angle of padeye load (deg)
    zlug : float
        Depth of the padeye (used for marker placement)
    z0 : float
        Depth of the mudline
    '''
    # Validate inputs
    if not drag_values or not depth_values:
        raise ValueError('drag_values and depth_values must be non-empty.')
    if len(drag_values) != len(depth_values):
        raise ValueError('drag_values and depth_values must be the same length.')
    
    fig, ax = plt.subplots(figsize=(8, 6))
        
    # Normalize values
    max_Su  = max((layer.get('Su_top', 0) for layer in layers if layer['soil_type'] == 'clay'), default=1)
    max_phi = max((layer.get('phi_top', 0) for layer in layers if layer['soil_type'] == 'sand'), default=1)

    seen_labels = set()
    for layer in layers:
        z_top = layer['top']
        z_bot = layer['bottom']
        soil = layer['soil_type']

        if soil == 'clay':
            Su = layer.get('Su_top', 25)
            gamma = layer.get('gamma_top', 10)
            color_fill = plt.cm.Oranges(Su/max_Su)
            label_soil = f'Su = {Su:.0f} kPa, γ = {gamma:.1f} kN/m³'
        elif soil == 'sand':
            phi = layer.get('phi_top', 30)
            gamma = layer.get('gamma_top', 10)
            color_fill = plt.cm.YlOrBr(phi / max_phi)
            label_soil = f'ϕ = {phi:.0f}°, γ = {gamma:.1f} kN/m³'
        else:
            color = 'gray'
            label = soil.capitalize()

        if label_soil not in seen_labels:
            ax.axhspan(-z_bot, -z_top, color=color_fill, alpha=0.4, label=label_soil)
            seen_labels.add(label_soil)
        else:
            ax.axhspan(-z_bot, -z_top, color=color_fill, alpha=0.4)
    
    scale = 2e6  # Arrow scaling factor for better visual readability
    
    # Plot the inverse catenary profile
    ax.plot(drag_values, depth_values, color='k', label='Mooring line')

    # Load arrows
    ax.arrow(0, -layers[0]['top'],
             Tm*np.cos(np.deg2rad(thetam))/scale, Tm*np.sin(np.deg2rad(thetam))/scale,
             head_width=0.25, head_length=0.5, color='r', label='Mudline load')

    ax.arrow(drag_values[-1], depth_values[-1],
             Ta*np.cos(np.deg2rad(thetaa))/scale, Ta*np.sin(np.deg2rad(thetaa))/scale,
             head_width=0.25, head_length=0.5, color='g', label='Lug load')
    
    ax.plot(0, -layers[0]['top'], 'ro', zorder=5)

    if zlug is not None:
        ax.plot(drag_values[-1], -zlug, 'go', label=f'Padeye (zlug = {zlug:.2f} m)')
    
    # Add mudline and padeye markers
    ax.axhline(-layers[0]['top'], color='b', linestyle='--', lw=1.5, label=f'Mudline')

    # Annotate loads
    ax.annotate(f"{Tm/1e6:.2f} MN", (Tm*np.cos(np.deg2rad(thetam))/scale, 
                                     -layers[0]['top']), color='r')
    ax.annotate(f"{Ta/1e6:.2f} MN", (drag_values[-1] + Ta*np.cos(np.deg2rad(thetaa))/scale,
                                     depth_values[-1]), color='g')

    # Deduplicate legend entries
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), loc='lower left')
    
    ax.set_xlabel('Drag distance (m)')
    ax.set_ylabel('Embedded depth (m)')
    ax.set_title('Inverse Catenary in Layered Soil')
    ax.grid(True)
    ax.set_ylim(min(zlug - 10, min(depth_values) - 5), max(5, max(depth_values) + 5))
    ax.legend(loc='lower right')
    plt.tight_layout()
    plt.show()

def plot_pycurve(pycurve_data):
    '''
    Plot multiple p–y curves from a mixed soil profile.

    Parameters
    ----------
    pycurve_data : list of tuples
        Each tuple must be (y_vals, p_vals, z_depth, soil_type)
    '''
    fig, ax = plt.subplots(figsize=(6, 5))

    for y, p, z, soil in pycurve_data:
        label = f'{soil.capitalize()} @ z = {z:.1f} m'
        ax.plot(y, p, label=label)

    ax.set_xlabel('Lateral displacement y (m)')
    ax.set_ylabel('Soil resistance p (N/m)')
    ax.set_title('p–y Curves at Various Depths')
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.show()
