
from anchor_map import Anchor
from famodel.anchors.anchors_famodel_map.capacity_plots_map import plot_load

# --- Define soil profile ---
profile_map = [
    {
        'name': 'CPT_T1',
        'x': 0.0, 'y': 0.0,
        'layers': [
            {'top':  0.0, 'bottom': 20.0, 'soil_type': 'clay', 'gamma_top': 8.0, 'gamma_bot': 8.5, 'Su_top':  50, 'Su_bot':  70},
            {'top': 20.0, 'bottom': 25.0, 'soil_type': 'clay', 'gamma_top': 8.5, 'gamma_bot': 8.5, 'Su_top':  80, 'Su_bot': 100},
            {'top': 25.0, 'bottom': 50.0, 'soil_type': 'clay', 'gamma_top': 8.5, 'gamma_bot': 9.0, 'Su_top': 125, 'Su_bot': 150}
        ]
    }
]

# Define the torpedo anchor
anchor = Anchor(
    dd = {
        'type': 'torpedo',
        'design': {
            'D1': 2.0,     # Wing diameter
            'D2': 1.5,     # Shaft diameter
            'L1': 11.0,    # Winged section length
            'L2': 5.0,     # Shaft section length
            'zlug': 20.0,  # Padeye depth
            'ballast': 10000
        }
    },
    r = [0.0, 0.0, 0.0]
)

# Assign applied loads
anchor.loads = {
    'Hm': 6.0e6,
    'Vm': 8.0e6
}

anchor.line_type = 'chain'
anchor.d = 0.16
anchor.w = 5000.0

# Assign soil profile
anchor.setSoilProfile(profile_map)

# Compute lug forces
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

# Compute anchor capacity
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
for key, value in anchor.capacity_results.items():
    print(f'{key}: {value:.2f}')


# --- Step 3: Optimize Anchor Geometry ---
anchor.getSizeAnchor(
    geom = [
        anchor.dd['design']['L1'], 
        anchor.dd['design']['D1']
        ],
    geomKeys = ['L1', 'D1'],
    geomBounds = [(7.0, 25.0), (2.5, 4.5)],
    loads = None,
    minfs = {'Ha': 1.0, 'Va': 1.0},
    lambdap_con = [2, 8],
    zlug_fix = True,
    plot = True
)

print('\nFinal Optimized Anchor:')
print('Design:', anchor.dd['design'])
print('Capacity Results:', anchor.capacity_results)

# --- Step 4: Visualize Anchor Geometry ---
# anchor.getCombinedPlot()
