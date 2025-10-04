
from famodel.anchors.anchor import Anchor

# --- Define soil profile ---
profile_map = [
    {
        'name': 'CPT_D1',
        'x': 0.0, 'y': 0.0,
        'layers': [
            {'top':  1.5, 'bottom':  5.0, 'soil_type': 'rock', 'UCS_top':  6.0, 'UCS_bot':  8.0, 'Em_top':  175, 'Em_bot':  290},
            {'top':  5.0, 'bottom':  9.0, 'soil_type': 'rock', 'UCS_top':  8.0, 'UCS_bot': 10.7, 'Em_top':  277, 'Em_bot':  297},
            {'top':  9.0, 'bottom': 30.0, 'soil_type': 'rock', 'UCS_top':  8.0, 'UCS_bot': 10.5, 'Em_top':  280, 'Em_bot': 305}
        ]
    }
]

# --- Create driven pile anchor ---
anchor = Anchor(
    dd = {
        'type': 'dandg',
        'design': {
            'L': 10.0,        # Embedded length
            'D': 2.85,        # Diameter
            'zlug': 1.0       # Padeye depth
        }
    },
    r = [0.0, 0.0, 0.0]
)

# Assign mooring loads
anchor.loads = {
    'Hm': 5.0e6,
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
    geomBounds = [(2.0, 30.0), (2.25, 5.0)],
    loads = None,
    lambdap_con = [4, 50],
    zlug_fix = True,
    safety_factor = {'SF_horizontal': 1, 'SF_vertical': 1},
    plot = True)