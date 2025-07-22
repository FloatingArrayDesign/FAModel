
import sys
sys.path.append(r'C:\Code\FAModel_anchors\famodel')

from project import Project
from anchors.anchor import Anchor

# Step 1: Initialize and load soil
proj = Project()
proj.loadSoil(
    filename='inputs/GulfOfMaine_soil_layered_100x100.txt',
    soil_mode='layered',
    profile_source='inputs/GulfOfMaine_soil_profiles.yaml')

for label, props in proj.soilProps.items():
    print(f"{label}: {props}")

# Step 2: Create and register an anchor at a known position in the grid
anchor = Anchor(
    dd = {'type': 'suction', 'design': {'D': 2.5, 'L': 15.0, 'zlug': 10.67}},
    r  = [54.0, -4450.0, 0.0])

# Step 3: Assign local soil profile from project (nearest neighbor lookup)
soil_id, soil_profile = proj.getSoilAtLocation(anchor.r[0], anchor.r[1])
anchor.soilProps = {soil_id: soil_profile}
anchor.setSoilProfile([{ 'name': soil_id, 'layers': soil_profile }])  # ensures `anchor.soil_profile` is set

# Step 4: Assign loads and line
anchor.loads = {'Hm': 2e6, 'Vm': 1.5e6}
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
anchor.getCostAnchor()
print(f'Material cost: {anchor.cost["Material cost"]:.2f} USD [2024]')

results = anchor.getSizeAnchor_gradient(
    geom=[anchor.dd['design']['L'], anchor.dd['design']['D']], 
    geomKeys= ['L','D'],
    geomBounds=[(12.0, 18.0), (1.5, 3.5)],
    safety_factor={'SF_combined': 1}, 
    zlug_fix=False, 
    lambdap_con=[4, 6], 
    step_size=0.2,
    tol=0.05,
    max_iter=30,
    verbose=True)

# Step 6: Report
print('\nFinal Optimized Anchor:')
print('Design:', anchor.dd['design'])
print('Capacity Results:', anchor.anchorCapacity)
anchor.getCostAnchor()
print(f'Material cost: {anchor.cost["Material cost"]:.2f} USD [2024]')


