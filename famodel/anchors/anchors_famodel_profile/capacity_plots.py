
import numpy as np
import matplotlib.pyplot as plt

def plot_pile(profile, soil_type, y, z, D, L, z0=None, zlug=None, hinge_location=None):
    '''Plots the deflected shape of the pile alongside the original vertical line.

    Parameters:
    ----------
    profile :
        
    soil_type : str
        
    y : np.array
        Lateral displacements (m)
    z : np.array
        Depths from pile head (m)
    D : float
        Pile diameter (m)
    L : float
        Pile length (m)
    z0 : float, optional
        Depth of the mudline from the pile head (m)
    zlug : float, optional
        Depth of the padeye from the pile head (m)
    hinge_location : int, optional
        Node index where a plastic hinge formed (if any)
    label : str
        Label for deflected shape line
    color : str
        Line color for deflected shape
    '''
    fig, ax = plt.subplots(figsize=(3, 5))
    
    lambdap = L / D

    # Adjust horizontal scale based on slenderness
    if lambdap >= 4:
        xmax = 5*D     # Slender (e.g., driven pile)
    elif lambdap <= 2:
        xmax = 2*D     # Stubby (e.g., drilled & grouted)
    else:
        xmax = 3*D     # Intermediate case

    # Mudline marker 
    if z0 is not None:
        ax.axhline(z0, color='blue', linestyle='--', label=f'Mudline (z0 = {z0:.2f} m)')
        #ax.fill_betweenx(z, -0.05 * D, 0.05 * D, where=(z >= z0), color='lightgray', alpha=0.3, label='Soil zone')

    # Padeye marker (on right wall of pile)
    if zlug is not None:
        ax.plot(D/2, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')

    # Draw pile as rectangle (from head to tip)
    pile = plt.Rectangle((-D/2, 0), D, L, edgecolor='k', facecolor='none', lw=2, label='Driven Pile')
    ax.add_patch(pile)
       
    seen_labels = set()
    # Plot soil layers as background fills
    for i in range(len(profile) - 1):
        z_top = profile[i][0]
        z_bot = profile[i+1][0]
        if soil_type == 'clay':
            Su = profile[i][1]
            color = plt.cm.Oranges(Su/np.max(profile[:, 1]))
            label = f'Su = {Su:.0f} kPa'
        elif soil_type == 'sand':
            phi = profile[i][1]
            gamma = profile[i][2]
            color = plt.cm.YlOrBr(phi/np.max(profile[:, 1]))
            label = f'ϕ = {phi:.0f}°, γ = {gamma:.1f} kN/m³'
        elif soil_type == 'rock':
            UCS = profile[i][1]
            Em = profile[i][2]
            color = plt.cm.Greys(UCS/np.max(profile[:, 1]))
            label = f'UCS = {UCS:.2f} MPa, Em = {Em:.1f} MPa'
            
        # Only assign label if not already used
        if label not in seen_labels:
            ax.axhspan(z_bot, z_top, color=color, alpha=0.4, label=label)
            seen_labels.add(label)
        else:
            ax.axhspan(z_bot, z_top, color=color, alpha=0.4)

    ax.set_xlabel('Lateral displacement (m)')
    ax.set_ylabel('Depth (m)')
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim([L + 5, -2])  # Downward is positive z
    ax.grid(ls='--')
    ax.legend()
    ax.set_title('Pile Deflection Profile')
    plt.tight_layout()
    plt.show()
    
    
def plot_suction(profile, soil_type, L, D, zlug=None, title='Suction Pile and Soil Layers', ax=None):
    '''Plot the soil profile and a suction pile geometry using array structure for profile.

    Parameters:
    ----------
        profile : list or ndarray of lists (z, Su or phi or UCS, gamma)
        soil_type : str
            'clay', 'sand', or 'rock'
        L : float
            Embedded length (m)
        D : float
            Pile diameter (m)
        zlug : float
            Padeye depth (m, referenced to pile head = 0)
        title : string
            Plot title
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    profile = np.array(profile)
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    xmax = 2*D

    # Split numerical values from profile
    z_vals = [float(row[0]) for row in profile]
    values = [float(row[1]) for row in profile]
    
    seen_labels = set()
    # Plot soil layers as background fills
    for i in range(len(z_vals) - 1):
        z_top = z_vals[i]
        z_bot = z_vals[i+1]
        val = values[i]

        if soil_type == 'clay':
            color = plt.cm.Oranges(val / max(values))
            label = f'Su = {val:.0f} kPa'
        elif soil_type == 'sand':
            color = plt.cm.YlOrBr(val / max(values))
            label = f'ϕ = {val:.0f}°'

        if label not in seen_labels:
            ax.axhspan(z_bot, z_top, color=color, alpha=0.4, label=label)
            seen_labels.add(label)
        else:
            ax.axhspan(z_bot, z_top, color=color, alpha=0.4)

    # Draw pile geometry
    x_left = -D/2; x_right = D/2
    z_top = 0; z_bot = L
    ax.plot([x_left,  x_left], [z_top, z_bot], color='k', lw=2.0, label='Suction Pile')
    ax.plot([x_right, x_right], [z_top, z_bot], color='k', lw=2.0)
    ax.plot([x_left, x_right], [z_top, z_top], color='k', lw=2.0)

    # Reference lines
    ax.axhline(z_vals[0], color='k', linestyle='--', lw=1.5, label='Mudline')
    ax.axhline(L, color='b', linestyle='--', lw=1.5, label='Pile Tip')

    # Padeye marker
    if zlug is not None and 0 <= zlug <= L:
        ax.plot(x_right, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')

    ax.set_xlabel('Horizontal extent (m)')
    ax.set_ylabel('Depth (m)')
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(L + 2*D, -D)
    ax.set_title(title)
    ax.grid()
    ax.legend()
    plt.tight_layout()
    plt.show()


def plot_torpedo(profile, soil_type, D1, D2, L1, L2, zlug, title='Torpedo Pile and Soil Layers'):
    '''Plot the soil layers and geometry of a torpedo pile using absolute depth for soil and pile head at z=0.

    Parameters:
    ----------
        profile : 
            list or array of (z, Su, 'clay')
        soil_type : str
            'clay'
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
    z1 = zlug + L1            # interface between L1 and L2
    z_tip = z1 + L2           # pile tip

    # Split numerical values from profile
    z_vals = [float(row[0]) for row in profile]
    values = [float(row[1]) for row in profile]

    seen_labels = set()
    # Plot soil layers as background fills
    for i in range(len(z_vals) - 1):
        z_top = z_vals[i]
        z_bot = z_vals[i+1]
        val = values[i]

        if soil_type == 'clay':
            color = plt.cm.Oranges(val/max(values))
            label = f'Su = {val:.0f} kPa'

        # Only assign label if not already used
        if label not in seen_labels:
            ax.axhspan(z_bot, z_top, color=color, alpha=0.4, label=label)
            seen_labels.add(label)
        else:
            ax.axhspan(z_bot, z_top, color=color, alpha=0.4)
        
    # Draw winged section (upper)
    ax.add_patch(plt.Rectangle((-D1/2, zlug), D1, L1, edgecolor='k', facecolor='none', lw=2, label='Winged Section'))

    # Draw shaft section (lower)
    ax.add_patch(plt.Rectangle((-D2/2, z1), D2, L2, edgecolor='k', facecolor='none', lw=2, label='Shaft Section'))
   
    # Reference lines
    ax.axhline(z_vals[0], color='k', linestyle='--', lw=1.0, label='Mudline')
    ax.axhline(    z_tip, color='b', linestyle='--', lw=1.0, label='Pile Tip')
    
    # Padeye marker
    if zlug is not None:
        ax.plot(0, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')


    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(max(z_vals) + 0.5*D1, min(zlug - 2*D1, z_vals[0] - 2))
    ax.set_xlabel('Horizontal extent (m)')
    ax.set_ylabel('Depth (m)')
    ax.set_title(title)
    ax.grid()
    ax.legend()
    plt.tight_layout()
    plt.show()
 
def plot_helical(profile, soil_type, D, L, d, zlug, n_helix=1, spacing=1.0, title='Helical Pile and Soil Layers'):
    '''Plot a helical pile in layered soil, with the pile starting at zlug and the helix near the pile tip.

    Parameters:
    ----------
        profile : list or array of (z, Su or phi, ...)
        soil_type : 
            'clay' or 'sand'
        D : float
            Helix diameter (m)
        L : float
            Pile length (m)
        d : float
            Shaft diameter (m)
        zlug : float
            Embedment depth of pile head (m)
        n_helix : int 
            Number of helices (-)
        spacing : float
            Vertical spacing between helices (m)
        title : str
            Plot title
    '''
    fig, ax = plt.subplots(figsize=(5, 5))
    
    xmax = 3*max(D, d)

    # Extract soil data
    z_vals = [float(row[0]) for row in profile]
    values = [float(row[1]) for row in profile]

    seen_labels = set()
    # Plot soil layers as background fills
    for i in range(len(z_vals) - 1):
        z_top = z_vals[i]
        z_bot = z_vals[i+1]
        val = values[i]
        
        if soil_type == 'clay':
            color = plt.cm.Oranges(val/max(values))
            label = f'Su = {val:.0f} kPa' 
        elif soil_type == 'sand':
            color = plt.cm.YlOrBr(val / max(values))
            label = f'ϕ = {val:.0f}°' 

        # Only assign label if not already used
        if label not in seen_labels:
            ax.axhspan(z_bot, z_top, color=color, alpha=0.4, label=label)
            seen_labels.add(label)
        else:
            ax.axhspan(z_bot, z_top, color=color, alpha=0.4)

    # Draw shaft
    ax.add_patch(plt.Rectangle((-d/2, 0), d, L, edgecolor='k', facecolor='none', lw=2, label='Shaft'))

    # Draw helices
    z_helix_base = L - D  # Base helix depth
    for i in range(n_helix):
        z_helix = z_helix_base - i*spacing
        if z_helix < zlug:
            break
        ax.plot([-D/2, D/2], [z_helix - d/2, z_helix + d/2], color='k', lw=2, label='Helix' if i == 0 else None)

    # Reference lines
    ax.axhline(z_vals[0], color='k', linestyle='--', lw=1.0, label='Mudline')
    ax.axhline(        L, color='b', linestyle='--', lw=1.0, label='Pile Tip')
    
    # Padeye marker
    if zlug is not None:
        ax.plot(d/2, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')

    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(max(z_vals), min(zlug - 0.5*d, z_vals[0] - 2))
    ax.set_xlabel('Horizontal extent (m)')
    ax.set_ylabel('Depth (m)')
    ax.set_title(title)
    ax.grid()
    ax.legend()
    plt.tight_layout()
    plt.show()

def plot_plate(profile, B, L, zlug, beta, title='Plate Anchor in Clay Profile'):
    '''Plot soil layers and an inclined plate anchor centered at zlug.

    Parameters:
    ----------
        profile : ndarray of shape (n, 3), with (depth, Su, gamma)
        B : float
            Plate width (m)
        L : float
            Plate length (m) 
        zlug : float
            Center embedment of the plate (m)
        beta : float
            Inclination angle of plate (deg)
        title : str
            Plot title
    '''
    fig, ax = plt.subplots(figsize=(5, 5))
    
    xmax = 3*B

    # Extract soil data
    layer_depths = profile[:, 0]
    layer_depths = np.append(layer_depths, [profile[0][-1]])    

    # Inclined plate geometry
    dx = (B/2)*np.cos(np.deg2rad(beta))
    dz = (B/2)*np.sin(np.deg2rad(beta))
    plate_x = [-dx, dx]
    plate_z = [zlug - dz, zlug + dz]

    seen_labels = set()
    # Plot soil layers as background fills
    for i in range(len(layer_depths) - 1):
        z_top = layer_depths[i]
        z_bot = layer_depths[i+1]
        Su = profile[i][1]
        color = plt.cm.Oranges(Su/np.max(profile[:, 1]))
        label = f'Su = {Su:.0f} kPa'
        
        # Only assign label if not already used
        if label not in seen_labels:
            ax.axhspan(z_bot, z_top, color=color, alpha=0.4, label=label)
            seen_labels.add(label)
        else:
            ax.axhspan(z_bot, z_top, color=color, alpha=0.4)

    # Plot inclined plate
    ax.plot(plate_x, plate_z, color='k', lw=1.5, label='Plate')

    # Reference lines
    ax.axhline(0, color='k', linestyle='--', lw=1.0, label='Mudline')
    
    # Padeye marker
    if zlug is not None:
        ax.plot(0, zlug, 'ko', label=f'Padeye (zlug = {zlug:.2f} m)')

    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(25, -1)  
    ax.set_xlabel("Horizontal extent (m)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(title)
    ax.legend(loc='lower right')
    ax.grid(True)
    plt.tight_layout()
    plt.show()
    
def plot_load(profile, soil_type, drag_values, depth_values, Tm, thetam, Ta, thetaa, zlug):

        fig, ax = plt.subplots(figsize=(12, 6))
        n = 2e6

        # Plot the inverse catenary profile
        ax.plot(drag_values, depth_values, color='b', label='Mooring line')

        # Add load vectors
        ax.arrow(0, -profile[0][0], Tm*np.cos(np.deg2rad(thetam))/n, Tm*np.sin(np.deg2rad(thetam))/n,
                 head_width=0.25, head_length=0.5, color='r', label='Mudline Load')
        ax.arrow(drag_values[-1], depth_values[-1], Ta*np.cos(thetaa)/n, Ta*np.sin(thetaa)/n,
                 head_width=0.25, head_length=0.5, color='g', label='Padeye Load')
        
        #ax.set_aspect('equal', adjustable='datalim')

        # Plot soil layers as background fills
        for i in range(len(profile) - 1):
            z_top = -profile[i][0]
            z_bot = -profile[i+1][0]
            if soil_type == 'clay':
                Su = profile[i][1]
                gamma = profile[i][2]
                color = plt.cm.Oranges(Su/np.max(profile[:, 1]))
                label = f'Su = {Su:.0f} kPa'
            elif soil_type == 'sand':
                phi = profile[i][1]
                gamma = profile[i][2]
                color = plt.cm.YlOrRd(phi/np.max(profile[:, 1]))
                label = f'ϕ = {phi:.0f}°, γ = {gamma:.1f} kN/m³'
            ax.axhspan(z_bot, z_top, color=color, alpha=0.4, label=label)

        # Deduplicate legend entries
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        ax.legend(unique.values(), unique.keys(), loc='lower left')

        ax.set_xlabel('Drag distance [m]')
        ax.set_ylabel('Embedded depth [m]')
        ax.set_title('Inverse Catenary in Layered Soil')
        ax.grid(True)
        
        ax.annotate(f"{Tm/1e6:.2f} MN", (Tm*np.cos(np.deg2rad(thetam))/n, -profile[0][0] + Tm*np.sin(np.deg2rad(thetam))/n), color='r')
        ax.annotate(f"{Ta/1e6:.2f} MN", (drag_values[-1] + Ta*np.cos(thetaa)/n, depth_values[-1] + Ta*np.sin(thetaa)/n), color='g')
