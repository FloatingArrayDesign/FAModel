
from famodel.anchors.anchor import Anchor
from famodel.anchors.anchors_famodel.support_plots import plot_suction

# --- Define soil profile ---
profile_map = [
    {
        'name': 'CPT_A1',
        'x': 0.0, 'y': 0.0,
        'layers': [
            {'top':  0.0, 'bottom': 12.0, 'soil_type': 'clay', 'gamma_top':  8.0, 'gamma_bot':  8.0, 'Su_top':  10, 'Su_bot':  20},
            {'top': 12.0, 'bottom': 22.0, 'soil_type': 'clay', 'gamma_top':  8.0, 'gamma_bot':  8.0, 'Su_top':  15, 'Su_bot':  25},
            {'top': 22.0, 'bottom': 30.0, 'soil_type': 'clay', 'gamma_top':  8.0, 'gamma_bot':  8.0, 'Su_top':  55, 'Su_bot':  70},
            {'top': 30.0, 'bottom': 40.0, 'soil_type': 'clay', 'gamma_top':  8.0, 'gamma_bot':  8.0, 'Su_top': 100, 'Su_bot': 100}
        ]
    },
    {
        'name': 'CPT_B1',
        'x': 500.0, 'y': 0.0,
        'layers': [
            {'top':  0.0, 'bottom':  5.0, 'soil_type': 'clay', 'gamma_top':  8.2, 'gamma_bot':  8.2, 'Su_top':  12, 'Su_bot':  22},
            {'top':  5.0, 'bottom': 15.0, 'soil_type': 'clay', 'gamma_top':  8.2, 'gamma_bot':  8.2, 'Su_top':  22, 'Su_bot':  22},
            {'top': 15.0, 'bottom': 30.0, 'soil_type': 'clay', 'gamma_top':  8.2, 'gamma_bot':  8.2, 'Su_top':  22, 'Su_bot':  50},
            {'top': 30.0, 'bottom': 40.0, 'soil_type': 'clay', 'gamma_top':  8.2, 'gamma_bot':  8.2, 'Su_top': 100, 'Su_bot': 100}
        ]
    },
    {
        'name': 'CPT_A2',
        'x': 0.0, 'y': 500.0,
        'layers': [
            {'top':  0.0, 'bottom':  8.0, 'soil_type': 'clay', 'gamma_top':  8.4, 'gamma_bot':  8.4, 'Su_top':  14, 'Su_bot':  14},
            {'top':  8.0, 'bottom': 18.0, 'soil_type': 'clay', 'gamma_top':  8.4, 'gamma_bot':  8.4, 'Su_top':  15, 'Su_bot':  45},
            {'top': 18.0, 'bottom': 30.0, 'soil_type': 'clay', 'gamma_top':  8.4, 'gamma_bot':  8.4, 'Su_top':  55, 'Su_bot':  55},
            {'top': 30.0, 'bottom': 40.0, 'soil_type': 'clay', 'gamma_top':  8.4, 'gamma_bot':  8.4, 'Su_top': 100, 'Su_bot': 100}
        ]
    },
    {
        'name': 'CPT_B2',
        'x': 500.0, 'y': 500.0,
        'layers': [
            {'top':  0.0, 'bottom': 15.0, 'soil_type': 'clay', 'gamma_top':  8.6, 'gamma_bot':  9.6, 'Su_top':  20, 'Su_bot':  26},
            {'top': 15.0, 'bottom': 25.0, 'soil_type': 'clay', 'gamma_top':  8.6, 'gamma_bot':  8.6, 'Su_top':  20, 'Su_bot':  40},
            {'top': 25.0, 'bottom': 30.0, 'soil_type': 'clay', 'gamma_top':  8.6, 'gamma_bot':  8.6, 'Su_top':  40, 'Su_bot':  40},
            {'top': 30.0, 'bottom': 40.0, 'soil_type': 'clay', 'gamma_top':  8.6, 'gamma_bot':  8.6, 'Su_top': 100, 'Su_bot': 100}
        ]
    }
]


anchor = Anchor(
    dd = {'type': 'suction', 'design': {'D': 3.5, 'L': 12.0, 'zlug': 8.67}},
    r = [250.0, 250.0, 000.0]
)

# --- Step 0: Create anchor based grid CPTs ---
anchor.interpolateSoilProfile(profile_map)

# --- Step 1: Plot suction pile and soil profile ---
# Access anchor geometrical properties
L = anchor.dd['design']['L']
D = anchor.dd['design']['D']
zlug = anchor.dd['design']['zlug']
# Access matched profile
layers = anchor.soil_profile
z0 = layers[0]['top']  

plot_suction(layers, L=L, D=D, z0=z0, zlug=zlug)

# Assign loads manually
anchor.loads = {
    'Hm': 3.68e6,    # Horizontal mudline load (N)
    'Vm': 0          # Vertical mudline load (N)
}

# Assign line properties manually
anchor.line_type = 'chain'
anchor.d = 0.16    # Chain diameter (m)
anchor.w = 5000.0  # Nominal submerged weight (N/m)


# --- Step 2: Compute Lug Forces ---
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

# --- Step 3: Compute Capacity ---
anchor.getCapacityAnchor(
    Hm = anchor.loads['Hm'],
    Vm = anchor.loads['Vm'],
    zlug = anchor.dd['design']['zlug'],
    line_type = anchor.line_type,
    d = anchor.d, w = anchor.w,
    mass_update=False,
    plot = True
)

print('\nCapacity Results:')
for key, value in anchor.anchorCapacity.items():
    print(f'{key}: {value:.2f}')
    
# --- Step 4: Compute Costs ---     
anchor.getCostAnchor()  

print(f"Mass: {anchor.anchorCapacity['Weight pile']/9.81:.2f} kg")
print(f"Material unit cost: {anchor.cost['unit_cost']:.2f} USD/kg")
print(f'Material cost: {anchor.cost["Material cost"]:.2f} USD [2024]')


# --- Step 5: Optimize Anchor Geometry ---
anchor.getSizeAnchor(
    geom = [anchor.dd['design']['L'], anchor.dd['design']['D']],
    geomKeys = ['L', 'D'],
    geomBounds = [(5.0, 15.0), (1.0, 4.0)],
    loads = None,
    lambdap_con = [3, 6],
    zlug_fix = False,
    safety_factor = {'SF_combined': 1},
    plot = True   
)

print('\nFinal Optimized Anchor:')
print('Design:', anchor.dd['design'])
print('Capacity Results:', anchor.anchorCapacity)

# # --- Step 6: Compute Costs ---     
anchor.getCostAnchor()  

print(f"Mass: {anchor.anchorCapacity['Weight pile']/9.81:.2f} kg")
print(f"Material unit cost: {anchor.cost['unit_cost']:.2f} USD/kg")
print(f'Material cost: {anchor.cost["Material cost"]:.2f} USD [2024]')

# --- Step 7: Visualize Anchor Geometry ---
# anchor.getCombinedPlot()
