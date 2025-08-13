
from famodel.anchors.anchor import Anchor

# --- Soil profile for helical pile in clay ---
profile_map = [
    {
        'name': 'CPT_H1',
        'x': 0.0, 'y': 0.0,
        'layers': [
            {'top': 1.0, 'bottom':  3.0, 'soil_type': 'clay', 'gamma_top':  8.0, 'gamma_bot':  9.0, 'Su_top':  60, 'Su_bot':  50},
            {'top': 3.0, 'bottom':  7.0, 'soil_type': 'clay', 'gamma_top': 15.0, 'gamma_bot': 25.0, 'Su_top': 100, 'Su_bot': 150},
            {'top': 7.0, 'bottom': 25.0, 'soil_type': 'clay', 'gamma_top': 25.0, 'gamma_bot': 50.0, 'Su_top': 200, 'Su_bot': 400}
        ]
    }
]

# --- Define helical anchor ---
anchor = Anchor(
    dd = {
        'type': 'helical',
        'design': {
            'D': 1.5,         # Helix diameter (m)
            'L': 12.0,        # Total length (m)
            'd': 0.5,         # Shaft diameter (m)
            'zlug': 3.0       # Padeye depth (m)
        }
    },
    r = [0.0, 0.0, 0.0]
)

# --- Assign mooring loads and properties ---
anchor.loads = {
    'Hm': 8e5,
    'Vm': 1e5}

anchor.line_type = 'chain'
anchor.d = 0.16
anchor.w = 5000.0

# --- Assign local soil ---
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
    plot = False)

print('\nCapacity Results:')
for key, val in anchor.anchorCapacity.items():
    print(f'{key}: {val:.2f}')

# --- Step 3: Optimize Anchor Geometry ---
anchor.getSizeAnchor(
    geom = [anchor.dd['design']['L'], anchor.dd['design']['D']],
    geomKeys = ['L', 'D'],
    geomBounds = [(6.0, 25.0), (0.5, 2.0)],
    loads = None,
    lambdap_con = [6, 15],
    zlug_fix = True,
    safety_factor = {'SF_horizontal': 1, 'SF_vertical': 1},   
    plot = False)
