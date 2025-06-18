
from famodel.anchors.anchor_profile import Anchor
from famodel.anchors.anchors_famodel_profile.capacity_plots import plot_load
import numpy as np  

# Define the soil profile
soil_profile = np.array([
    [ 1.0,  10, 8.0],
    [ 2.0,  25, 8.5],
    [ 8.0,  50, 9.0],
    [16.0, 100, 9.5],
    [25.0, 100, 9.5]
])

# Create Anchor object
anchor = Anchor(
    dd={
        'type': 'suction',
        'design': {'D': 2.5, 'L': 10.0, 'zlug': 6.0, 'soil_type': 'clay'},
        'soil_properties': {'clay': soil_profile}
    },
    ms=None,
    r=[0.0, 0.0, 0.0],
    aNum=0,
    id='A1',
    g=9.81,
    rho=1025
)

# Assign loads manually
anchor.loads = {
    'Hm': 3e6,     # Horizontal mudline load (N)
    'Vm': 1e6      # Vertical mudline load (N)
}

# Also assign mooring line properties manually
anchor.line_type = 'chain'
anchor.d = 0.16    # Chain diameter (m)
anchor.w = 5000.0  # Nominal submerged weight (N/m)

# --- Step 1: Compute Lug Forces ---
Ha, Va = anchor.getLugForces(
    ground_conds=anchor.dd['soil_properties'],
    Hm=anchor.loads['Hm'],
    Vm=anchor.loads['Vm'],
    thetam=np.degrees(np.arctan2(anchor.loads['Vm'], anchor.loads['Hm'])),
    zlug=anchor.dd['design']['zlug'],
    line_type=anchor.line_type,
    d=anchor.d,
    w=anchor.w,
    plot=True  
)

# Print Lug Forces
print('\nLug Forces Computed:')
print(f'Ha = {Ha:.2f} N')
print(f'Va = {Va:.2f} N')

# --- Step 2: Compute Anchor Capacity ---
anchor.getCapacityAnchor(
    ground_conds=anchor.dd['soil_properties'],
    Hm=anchor.loads['Hm'],
    Vm=anchor.loads['Vm'],
    thetam=np.degrees(np.arctan2(anchor.loads['Vm'], anchor.loads['Hm'])),
    zlug=anchor.dd['design']['zlug'],
    line_type=anchor.line_type,
    d=anchor.d,
    w=anchor.w,
    plot=True  
)

# Print Capacity Results
print('\nCapacity Results:')
for key, value in anchor.capacity_results.items():
    print(f'{key}: {value:.2f}')
    
# --- Step 3: Optimize Anchor Geometry ---
anchor.getSizeSuction(
    geom=[anchor.dd['design']['L'], anchor.dd['design']['D']],
    geomKeys=['L', 'D'],
    geomBounds=[(5.0, 15.0), (2.0, 6.0)],
    loads=None,
    minfs={'Ha': 1.0, 'Va': 1.0},
    lambdap_con=[3, 6],
    zlug_fix=False,
    plot=True
)

print('\nFinal Optimized Anchor:')
print('Design:', anchor.dd['design'])
print('Capacity Results:', anchor.capacity_results)

# --- Step 4: Visualize Anchor Geometry ---
anchor.getCombinedPlot()