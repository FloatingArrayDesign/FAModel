
from famodel.anchors.anchor import Anchor

# --- Define soil profile ---
profile_map = [
    {
        'name': 'CPT_D1',
        'x': 0.0, 'y': 0.0,
        'layers': [
            {'top':  1.5, 'bottom':  6.0, 'soil_type': 'rock', 'UCS_top':  5.0, 'UCS_bot':  5.0, 'Em_top':  7, 'Em_bot': 7},
            {'top':  6.0, 'bottom': 15.0, 'soil_type': 'rock', 'UCS_top':  6.0, 'UCS_bot':  6.7, 'Em_top':  7, 'Em_bot': 7},
            {'top': 15.0, 'bottom': 35.0, 'soil_type': 'rock', 'UCS_top': 10.0, 'UCS_bot': 10.5, 'Em_top':  7, 'Em_bot': 7}
        ]
    }
]

# --- Create driven pile anchor ---
anchor = Anchor(
    dd = {
        'type': 'driven',
        'design': {
            'L': 15.0,        # Embedded length
            'D': 1.85,        # Diameter
            'zlug': 1.5       # Padeye depth
        }
    },
    r = [0.0, 0.0, 0.0]
)

# Assign mooring loads
anchor.loads = {
    'Hm': 2.5e6,
    'Vm': 2.5e5}

anchor.line_type = 'chain'
anchor.d = 0.16
anchor.w = 5000.0

# Assign local soil
anchor.setSoilProfile(profile_map)

# --- Step 1: Capacity ---
anchor.getCapacityAnchor(
    Hm = anchor.loads['Hm'],
    Vm = anchor.loads['Vm'],
    zlug = anchor.dd['design']['zlug'],
    line_type = anchor.line_type,
    d = anchor.d,
    w = anchor.w,
    plot = True)

print('\nCapacity Results:')
for key, val in anchor.anchorCapacity.items():
    print(f'{key}: {val:.2f}')

# --- Step 2: Optimize Anchor Geometry ---
anchor.getSizeAnchor(
    geom = [anchor.dd['design']['L'], anchor.dd['design']['D']],
    geomKeys = ['L', 'D'],
    geomBounds = [(2.0, 70.0), (0.25, 3.0)],
    loads = None,
    lambdap_con = [4, 50],
    zlug_fix = True,
    safety_factor = {'SF_horizontal': 1, 'SF_vertical': 1},
    plot = True)