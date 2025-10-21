
from famodel.anchors.anchor import Anchor

# --- Define soil profile ---
profile_map = [
    {
        'name': 'CPT_D1',
        'x': 0.0, 'y': 0.0,
        'layers': [
            {'top':  1.5, 'bottom':  6.0, 'soil_type': 'clay', 'gamma_top':  9.0, 'gamma_bot': 10.0, 'Su_top':  25, 'Su_bot': 200},
            {'top':  6.0, 'bottom': 15.0, 'soil_type': 'sand', 'gamma_top': 10.0, 'gamma_bot': 10.0, 'phi_top': 28, 'phi_bot': 32, 'Dr_top':  80, 'Dr_bot': 85},
            {'top': 15.0, 'bottom': 35.0, 'soil_type': 'clay', 'gamma_top': 10.0, 'gamma_bot': 10.5, 'Su_top': 100, 'Su_bot': 100}
        ]
    }
]

# --- Create driven pile anchor ---
anchor = Anchor(
    dd = {
        'type': 'driven',
        'design': {
            'L': 25.0,        # Embedded length
            'D': 2.00,        # Diameter
            'zlug': 3.0       # Padeye depth
        }
    },
    r = [0.0, 0.0, 0.0]
)

# Assign mooring loads
anchor.loads = {
    'Hm': 7.0e5,
    'Vm': 2.5e5}

anchor.line_type = 'chain'
anchor.d = 0.16
anchor.w = 5000.0

# Assign local soil
anchor.setSoilProfile(profile_map)

# --- Step 1: Lug Forces ---
layers, Ha, Va = anchor.getLugForces(
    Hm = anchor.loads['Hm'],
    Vm = anchor.loads['Vm'],
    zlug = anchor.dd['design']['zlug'],
    line_type = anchor.line_type,
    d = anchor.d,
    w = anchor.w,
    plot = True)

print('\nLug Forces Computed:')
print(f'Ha = {Ha:.2f} N')
print(f'Va = {Va:.2f} N')

# --- Step 2: Capacity ---
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

# --- Step 3: Optimize Anchor Geometry ---
anchor.getSizeAnchor(
    geom = [anchor.dd['design']['L'], anchor.dd['design']['D']],
    geomKeys = ['L', 'D'],
    geomBounds = [(2.0, 50.0), (0.25, 3.5)],
    loads = None,
    lambdap_con = [3, 50],
    zlug_fix = True,
    safety_factor = {'SF_horizontal': 1, 'SF_vertical': 1},
    plot = True)