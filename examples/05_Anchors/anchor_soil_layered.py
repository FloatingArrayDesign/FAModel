

from famodel import Project
from famodel.anchors.anchor import Anchor

# Step 1: Initialize and load soil
proj = Project()
proj.loadSoil(
    filename='inputs/GulfOfMaine_soil_layered_100x100.txt',
    soil_mode='layered',
    profile_source='inputs/GulfOfMaine_soil_profiles.yaml')

# Step 2: Create and register an anchor at a known position in the grid
anchor = Anchor(
    dd = {'type': 'suction', 'design': {'D': 3.5, 'L': 12.0, 'zlug': 9.67}},
    r  = [54.0, -4450.0, 0.0])

# Step 3: Assign local soil profile from project (nearest neighbor lookup)
proj.getSoilAtLocation(anchor.r[0], anchor.r[1])
anchor.setSoilProfile(proj.profile_map)  

# Step 4: Assign loads and line
anchor.loads = {'Hm': 1e6, 'Vm': 1.5e6}
anchor.line_type = 'chain'
anchor.d = 0.16
anchor.w = 5000.0

# Step 5: Run capacity check and optimization
anchor.getLugForces(
    Hm=anchor.loads['Hm'], Vm=anchor.loads['Vm'], 
    zlug = anchor.dd['design']['zlug'], 
    d=anchor.d, w=anchor.w,  
    plot=True)

anchor.getCapacityAnchor(
    Hm=anchor.loads['Hm'], Vm=anchor.loads['Vm'], 
    zlug = anchor.dd['design']['zlug'], 
    line_type=anchor.line_type, d=anchor.d, w=anchor.w,   
    mass_update=True,
    plot=True)
anchor.getCost()
print(f'Material cost: {anchor.cost["Material cost"]:.2f} USD [2024]')

results = anchor.getSizeAnchor(
    geom = [anchor.dd['design']['L'], anchor.dd['design']['D']],
    geomKeys = ['L', 'D'],
    geomBounds = [(8.0, 15.0), (2.0, 4.0)],
    loads = None,
    lambdap_con = [3, 6],
    zlug_fix = False,
    safety_factor = {'SF_combined': 2},
    plot = True)

# Step 6: Report
print('\nFinal Optimized Anchor:')
print('Design:', anchor.dd['design'])
print('Capacity Results:', anchor.anchorCapacity)
anchor.getCost()
print(f'Material cost: {anchor.cost["Material cost"]:.2f} USD [2024]')


