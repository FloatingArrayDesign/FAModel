
from anchor_map import Anchor

# --- Define soil profile ---
profile_map = [
    {
        'name': 'CPT_H1',
        'x': 0.0, 'y': 0.0,
        'layers': [
            {'top':  1.0, 'bottom': 10.0, 'soil_type': 'sand', 'gamma_top': 10.0, 'gamma_bot': 11.0, 'phi_top': 30, 'phi_bot': 32, 'Dr_top': 60, 'Dr_bot': 60},
            {'top': 10.0, 'bottom': 20.0, 'soil_type': 'sand', 'gamma_top': 11.0, 'gamma_bot': 11.5, 'phi_top': 36, 'phi_bot': 38, 'Dr_top': 60, 'Dr_bot': 80}
        ]
    }
]

# --- Create helical anchor object ---
anchor = Anchor(
    dd = {
        'type': 'helical',
        'design': {
            'D': 1.7,         # Helix diameter (m)
            'L': 12.0,        # Depth (m)
            'd': 0.3,         # Shaft diameter (m)
            'zlug': 4.0       # Padeye depth (m)
        }
    },
    r = [0.0, 0.0, 0.0]
)

# Assign loads and mooring info
anchor.loads = {
    'Hm': 8.8e6,
    'Vm': 1.2e6
}
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
    plot = True
)

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
    plot = True
)

print('\nCapacity Results:')
for key, val in anchor.capacity_results.items():
    print(f'{key}: {val:.2f}')
